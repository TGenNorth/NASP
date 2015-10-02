package command

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"
)

var cmdFrankenfasta = &Command{
	UsageLine: "frankenfasta [mummer deltas]",
	Short:     "1-to-1 position align sample with a reference fasta",
	Long: `
Frankenfasta creates a fasta 1-to-1 position aligned with a reference fasta
by combining a sample fasta with a nucmer alignment delta file.

How it works:

To start with, a list of reference fasta contigs and their lengths are read from
the delta file. This is used to create a fasta where all positions are masked
with 'X'.

	>Contig Description
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

The delta alignments are used to copy positions from the sample fasta into their
corresponding positions in the reference fasta.

	>Contig Description
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
	XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXAAAAAAAXAAXXAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
	AAAAAGTTACTTTCATTAATAGAGCAAAATTTATTAATTATACTTTTTACACAAAAAACAAAGAGAGGATCGACCTTTTC
	TATTCATTCTTTAAACTACAGATTCTTTATAAAAAAATTTAACGACGTTGAGAGTAGACATGAGCGAAAGCACGGGTAGT
	AGCACGCTTCAACTTGGATTGAGCGGCCTTGAGGTCATCAGCTTGAGCATCGTCAGCGGTGAGTTCCATCAATTCCATGG
	CTTCAGCGAGAGCAGTTTCAATGGCTTCACGGTCACCACGACGGAGCTTCATGTTGACGTTAGGGTCAGAGATGTTGTTC
	TCAACTTGGTAGACGTAAGATTCAAGATCTTGCTTGGCTTGTACAACGGCTTCACGTTGCTTGTCAGATTCAGCATTCTT
	TTCAGCATCTTGAACCATACGTTCGATATCAGCAGCAGAGAGACGAGTAGAGTTGGAGATGGTGACATCAGCCTTGCGAC
	CAGTGGTCTTGTCTTGAGCAGTGACCTTCAAGATACCATTGGCATCCAAGTCAAAGGTACAGAGAAGTTCAGGGGTACCA
	CGGGGAGCAGGGATAAGACCAGAGAGAGAGAACTCACCGAGAAGGTTGTTTTCCTTGGTCATCAAACGTTCACCTTCGTA
	GACAGGGAAAGTGACAGTGGTTXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

'X' indicates the absence of calls. Either no sequence aligned to the position,
or it was a deletion.

'N' indicates uncertain calls. Two alignments with different calls may have mapped
to the same position.

Pre-release versions of NASP used a '.' to differentiate deletions from
	`,
}

var frankenfastaSample string // frankenfasta --sample flag

func init() {
	cmdFrankenfasta.Run = runFrankenfasta

	cmdFrankenfasta.Flag.StringVar(&frankenfastaSample, "sample", "", "")
}

func runFrankenfasta(cmd *Command, args []string) {
	for _, delta := range args {
		if err := frankenfasta(delta); err != nil {
			log.Println(err)
			cmd.Usage()
		}
	}
}

// TODO: document it will use the fasta named in the delta file
/*
Unimplemented: this template was part of an idea to output job scripts / job arrays
so the user could easily see / modify the commands that were run to produce the
nasp output rather than interpretting the run log.
const (
	convertExternalGenomeTemplate = `
external_fasta=${file[$index]}
nucmer --prefix ${sample_name} ${reference_fasta} ${external_fasta}
delta-filter -q -r -o 100 ${sample_name}.delta |
nasp frankenfasta - > ${sample_name}.frankenfasta
`
)
*/

var (
	// TODO: Report line number and content
	ErrUnexpectedFormat = errors.New("frankenfasta: parse error")
)

// TODO: The mummer package implements a delta parser. Delete this duplicate
// implementation and use the mummer package delta parser instead.
type header struct {
	referenceSequence string
	querySequence     string
	referenceLength   int64
	queryLength       int64
}

func (h *header) unmarshal(fields []string) (err error) {
	h.referenceSequence = fields[0][1:]
	h.querySequence = fields[1]
	if h.queryLength, err = strconv.ParseInt(fields[3], 10, 0); err != nil {
		return err
	}
	if h.referenceLength, err = strconv.ParseInt(fields[2], 10, 0); err != nil {
		return err
	}

	return nil
}

/*
Following this sequence header is the alignment data. Each alignment following
also has a header that describes the coordinates of the alignment and some
error information. These coordinates are inclusive and reference the forward
strand of the DNA sequence, regardless of the alignment type
(DNA or amino acid). Thus, if the start coordinate is greater than the end
coordinate, the alignment is on the reverse strand. The four coordinates are
the start and end in the reference and the start and end in the query
respectively. The three digits following the location coordinates are the
number of errors (non-identities + indels), similarity errors (non-positive
match scores), and stop codons (does not apply to DNA alignments, will be "0").
An example header might look like:
	2631 3401 2464 3234 15 15 2
*/
type alignment struct {
	referenceStart   int64
	referenceEnd     int64
	queryStart       int64
	queryEnd         int64
	numberOfErrors   int64
	similarityErrors int64
	stopCodons       int64
	distances        []int64
}

type deltaRecord struct {
	header     header
	alignments []alignment
}

func (d *deltaRecord) appendAlignment(fields []string) (err error) {
	var a alignment

	if a.referenceStart, err = strconv.ParseInt(fields[0], 10, 0); err != nil {
		return err
	}
	if a.referenceEnd, err = strconv.ParseInt(fields[1], 10, 0); err != nil {
		return err
	}
	if a.queryStart, err = strconv.ParseInt(fields[2], 10, 0); err != nil {
		return err
	}
	if a.queryEnd, err = strconv.ParseInt(fields[3], 10, 0); err != nil {
		return err
	}
	if a.numberOfErrors, err = strconv.ParseInt(fields[4], 10, 0); err != nil {
		return err
	}
	if a.similarityErrors, err = strconv.ParseInt(fields[5], 10, 0); err != nil {
		return err
	}
	if a.stopCodons, err = strconv.ParseInt(fields[6], 10, 0); err != nil {
		return err
	}

	d.alignments = append(d.alignments, a)

	return nil
}

func (d *deltaRecord) appendDistance(field string) (err error) {
	distance, err := strconv.ParseInt(field, 10, 0)
	if err != nil {
		return err
	}

	lastAlignment := len(d.alignments) - 1
	d.alignments[lastAlignment].distances = append(d.alignments[lastAlignment].distances, distance)

	return nil
}

type delta struct {
	dataType          string
	referenceFilepath string
	queryFilepath     string
	records           []*deltaRecord
}

/*
Below is an example of what a delta file might look like:
/home/username/reference.fasta /home/username/query.fasta
PROMER
>tagA1 tagB1 3000000 2000000
1667803 1667078 1641506 1640769 14 7 2
-145
-3
-1
-40
0
1667804 1667079 1641507 1640770 10 5 3
-146
-1
-1
-34
0
>tagA2 tagB4 4000 3000
2631 3401 2464 3234 4 0 0
0
2608 3402 2456 3235 10 5 0
7
1
1
1
1
0
(output continues ...)
*/
// FIXME: To fully implement the io.ReadFrom interface, it must return the
// number of bytes read from the io.Reader. Currently it returns nonsense.
func (d *delta) ReadFrom(r io.Reader) (n int64, err error) {
	var record *deltaRecord

	scanner := bufio.NewScanner(r)

	// The first line is the reference and query fasta filepaths separated by a space
	if !scanner.Scan() {
		return n, ErrUnexpectedFormat
	}

	// TODO: verify fields is composed of the reference and query filepaths
	fields := strings.Fields(scanner.Text())
	d.queryFilepath = fields[1]

	// The second line specifies the alignment data type
	if !scanner.Scan() {
		return n, ErrUnexpectedFormat
	}

	// TODO: Add PROMER?
	if d.dataType = scanner.Text(); d.dataType != "NUCMER" {
		return n, ErrUnexpectedFormat
	}

	// Parse the delta records
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())
		switch {
		default:
			return n, ErrUnexpectedFormat
		case strings.HasPrefix(fields[0], ">"):
			// Header line format: >reference_contig query_contig reference_length query_length
			if record != nil {
				// Append the previous record before starting a new one
				d.records = append(d.records, record)
			}
			record = &deltaRecord{}
			if record.header.unmarshal(fields); err != nil {
				return n, err
			}
		case len(fields) == 7:
			// Alignment line format: reference_start reference_end query_start query_end
			if err := record.appendAlignment(fields); err != nil {
				return n, err
			}
		case len(fields) == 1:
			// Distance is a cumulative offset to the next insert/delete
			if err := record.appendDistance(fields[0]); err != nil {
				return n, err
			}
		}
	}

	d.records = append(d.records, record)

	return n, nil
}

type fasta map[string][]byte

// Write contigs in fasta format to io.Writer wrapping sequences at 80 characters
func (f fasta) WriteTo(w io.Writer) {
	var keys []string

	for k := range f {
		keys = append(keys, k)
	}

	sort.Strings(keys)

	for _, k := range keys {
		sequence := f[k]

		fmt.Fprintf(w, ">franken::%s\n", k)
		for i := 80; i < len(sequence); i += 80 {
			fmt.Fprintf(w, "%s\n", sequence[i-80:i])
		}
		if remainder := len(sequence) % 80; remainder >= 0 {
			if remainder == 0 {
				remainder = 80
			}
			fmt.Fprintf(w, "%s\n", sequence[len(sequence)-remainder:])
		}
	}
}

func (f fasta) Unmarshal(r io.Reader, delta delta) {
	var isReferenceAligned bool
	var currentContig string
	var sequence []byte

	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		line := scanner.Bytes()
		if len(line) > 0 && line[0] == '>' {
			// Add the contig sequence before starting a new one
			if isReferenceAligned {
				f[currentContig] = sequence
			}
			// currentContig is the first word in a contig description.
			// Everything after the first space is stripped.
			// Given the following fasta contig description:
			// >gi|653474994|gb|ATWT01000581.1| Yersinia pestis EBD10-058 contig000581, whole genome shotgun sequence
			// currentContig = gi|653474994|gb|ATWT01000581.1|
			currentContig = string(bytes.SplitN(line[1:], []byte(" "), 2)[0])
			isReferenceAligned = false
			for i := range delta.records {
				if delta.records[i].header.querySequence == currentContig {
					isReferenceAligned = true
					sequence = make([]byte, 0, delta.records[i].header.queryLength)
					break
				}
			}
		} else if isReferenceAligned {
			sequence = sequence[0 : len(sequence)+len(line)]
			copy(sequence[len(sequence)-len(line):], line)
		}
		// else ignore contigs that did not align to the reference
	}
	if _, ok := f[currentContig]; !ok && isReferenceAligned {
		f[currentContig] = sequence
	}
}

func frankenfasta(deltaFilename string) error {
	var delta delta

	deltaFile, err := os.Open(deltaFilename)
	if err != nil {
		return err
	}
	defer deltaFile.Close()

	if _, err := delta.ReadFrom(deltaFile); err != nil {
		return err
	}

	fastaFile, err := os.Open(delta.queryFilepath)
	if err != nil {
		return err
	}
	defer fastaFile.Close()

	external := make(fasta)
	//log.Printf("Reading sample sequences from %s\n", delta.queryFilepath)
	external.Unmarshal(fastaFile, delta)

	franken := make(fasta)

	// Copy the aligned sample positions to their corresponding reference positions.
	// - Write a '.' when the position is a deletion
	// - Skip insertions
	// - Positions without an alignment will remain an 'X' indicating no data
	// The end result of a frankenfasta is a 1-to-1 position alignment of the sample to the reference
	for i := range delta.records {
		var frankenSequence []byte
		var isAlignmentOverlap bool

		if frankenSequence, isAlignmentOverlap = franken[delta.records[i].header.referenceSequence]; !isAlignmentOverlap {
			// Initialize frankenSequence with X
			frankenSequence = make([]byte, delta.records[i].header.referenceLength)
			for p := 0; p < len(frankenSequence); p++ {
				frankenSequence[p] = 'X'
			}
		}

		sequence, ok := external[delta.records[i].header.querySequence]
		if !ok {
			// TODO: replace with error type
			return fmt.Errorf("nasp.frankenfasta: the query sequence '%s' in %s was not found in %s\n", delta.records[i].header.querySequence, deltaFilename, delta.queryFilepath)
		}
		//log.Println(delta.records[i].header.querySequence)

		alignments := delta.records[i].alignments
		for _, alignment := range alignments {
			// Position 1-based indexing to array 0-based indexing
			refStart := alignment.referenceStart - 1
			qStart := alignment.queryStart - 1
			isQueryReversed := alignment.queryEnd < alignment.queryStart

			if isQueryReversed {
				qStart++
			}

			//log.Printf("%#v\n", alignment)

			//for dIdx, d := range alignment.distances {
			for _, d := range alignment.distances {
				var refEnd, qEnd int64
				var absDist int64

				// Distance is a cumulative offset to the next insert/delete.
				// The negative sign differentiates deletes from inserts.
				// 1 and -1 are short-circuit paths for a sequence of inserts/deletes.
				switch {
				case d == -1:
					// Skip the deleted position
					if isQueryReversed {
						qStart--
					} else {
						qStart++
					}
					//log.Println(dIdx, "d=", d, "refStart=", refStart, "qStart=", qStart, "isQueryReversed=", isQueryReversed)
					continue
				case d == 1:
					frankenSequence[refStart] = '.'
					//log.Println(dIdx, "d=", d, "refStart=", refStart, "qStart=", qStart, "isQueryReversed=", isQueryReversed)
					refStart++
					continue
				case d == 0:
					absDist = alignment.referenceEnd - refStart
				case d > 0:
					absDist = d - 1
				case d < 0:
					absDist = -d - 1
				}

				if absDist < 1 {
					// TODO: clarify when this condition is met
					// Occurs when the previous offset included the last position
					// and the 0 merely marks the end of the alignment
					break
				}

				refEnd = refStart + absDist

				if isQueryReversed {
					// Positions are read from end to beginning and reverse complemented
					qEnd = qStart - absDist
					//log.Println(dIdx, "d=", d, "absDist=", absDist, "refLen=", len(frankenSequence), "qLen=", len(sequence), "refStart=", refStart, "refEnd=", refEnd, "rDiff=", refEnd-refStart, "qStart=", qStart, "qEnd=", qEnd, "qDiff=", qStart-qEnd, "isQueryReversed=", isQueryReversed)
					segment := reverseComplement(sequence[qEnd:qStart])
					if isAlignmentOverlap {
						mergeMarkingConflictsWithN(frankenSequence[refStart:refEnd], segment)
					} else {
						copy(frankenSequence[refStart:refEnd], segment)
					}
				} else {
					// Positions are read from beginning to end
					qEnd = qStart + absDist
					//log.Println(dIdx, "d=", d, "absDist", absDist, "refLen=", len(frankenSequence), "qLen=", len(sequence), "refStart=", refStart, "refEnd=", refEnd, "rDiff=", refEnd-refStart, "qStart=", qStart, "qEnd=", qEnd, "qDiff=", qEnd-qStart, "isQueryReversed=", isQueryReversed)
					if isAlignmentOverlap {
						mergeMarkingConflictsWithN(frankenSequence[refStart:refEnd], sequence[qStart:qEnd])
					} else {
						copy(frankenSequence[refStart:refEnd], sequence[qStart:qEnd])
					}
				}

				switch {
				// Do nothing when the distance is 0. The positions in this region all mapped 1-to-1.
				//case d == 0:
				case d > 0:
					// Fill sample deletes with '.'
					frankenSequence[refEnd] = '.'
					refEnd++
				case d < 0:
					// Skip sample inserts
					if isQueryReversed {
						qEnd--
					} else {
						qEnd++
					}
				}

				refStart = refEnd
				qStart = qEnd

				//log.Println("end refStart=", refStart, "qStart=", qStart)
			}
		}

		franken[delta.records[i].header.referenceSequence] = frankenSequence
	}

	franken.WriteTo(os.Stdout)

	return nil
}

func reverseComplement(segment []byte) []byte {
	complement := make([]byte, len(segment))
	copy(complement, segment)
	l := len(complement) - 1
	h := l/2 + 1
	for i := 0; i <= l; i++ {
		// Swap positions until the halfway point to reverse the sequence
		if i < h {
			complement[i], complement[l-i] = complement[l-i], complement[i]
		}
		// Capitalize
		if 'a' <= complement[i] && complement[i] <= 'z' {
			complement[i] -= 'a' - 'A'
		}
		switch complement[i] {
		case 'A':
			complement[i] = 'T'
		case 'B':
			complement[i] = 'V'
		case 'C':
			complement[i] = 'G'
		case 'D':
			complement[i] = 'H'
		case 'G':
			complement[i] = 'C'
		case 'H':
			complement[i] = 'D'
		case 'M':
			complement[i] = 'K'
		case 'N':
			complement[i] = 'N'
		case 'R':
			complement[i] = 'Y'
		case 'S':
			complement[i] = 'S'
		case 'T':
			complement[i] = 'A'
		case 'U':
			complement[i] = 'A'
		case 'V':
			complement[i] = 'B'
		case 'W':
			complement[i] = 'W'
		case 'X':
			complement[i] = 'X'
		case 'Y':
			complement[i] = 'R'
		}
	}
	return complement
}

// mergeMarkingConflictsWithN merges 2 equal length DNA sequences.
// If two alignments overlap and they do not match perfectly, replace the
// conflicting calls with N.
func mergeMarkingConflictsWithN(dst, src []byte) {
	for i := range src {
		switch {
		// default: do nothing; both sequences match
		case dst[i] == 'X':
			// Fill unset/undefined positions with src calls
			dst[i] = src[i]
		case dst[i] != src[i]:
			// Replace conflicting calls with 'N'
			dst[i] = 'N'
		}
	}
}
