package command

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"sort"
	"strconv"
	"strings"

	"github.com/TGenNorth/nasp/mummer"
)

var cmdDuplicates = &Command{
	UsageLine: "duplicates --type [range|fasta] reference.delta > duplicates.txt",
	Short:     "identify duplicate regions in the reference",
	Long: `
Duplicates identifies duplicate regions in the reference genome contigs.

reference.delta is a delta encoded file from the nucmer aligner. It is created
by aligning a reference assembly against itself to identify duplicate regions.
	nucmer --maxmatch --nosimplify reference.fasta reference.fasta
The result is passed through the delta-filter pairing each position with its
best match.
	delta-filter -q -r out.delta > reference.delta

The --type flag sets the output type (default: fasta)

If --type is fasta, duplicates uses the delta to create a fasta-style output
where the nucleotide sequences are replaced with a binary sequence such that
'0' indicates a unique position and '1' indicates a duplicates position:

	>Contig Description
	00000000000000000000000000000000000000000000000000000000000000000000000000000000
	00000000000000000000000000000000000000000000000000000000000000000000000000000000
	00000000000000000000000000000000000000000000000000000000000000000000000000000000
	00000000000000000000000000000000000000000000000000000000000000000000000000000000
	00000000000000000000000000000000000000000000000000000000000000111111111111111111
	11111111111111111111111111111111111111111111111111111111111111111111111111111111
	11111111111111111111111111100

If --type is range, for each duplicate region duplicates outputs a line containing
the contig name, start position, and end position.

	Contig Description 500 1000
`,
}

var typeFlag string

func init() {
	cmdDuplicates.Run = runDuplicates
	cmdDuplicates.Flag.StringVar(&typeFlag, "type", "fasta", "")
}

func runDuplicates(cmd *Command, args []string) {
	if len(args) != 1 {
		log.Println("duplicates requires a delta file")
		cmdDuplicates.Usage()
	}

	delta, err := os.Open(args[0])
	if err != nil {
		log.Fatal(err)
	}
	defer delta.Close()

	br := bufio.NewReader(delta)
	bw := bufio.NewWriter(os.Stdout)
	defer bw.Flush()

	switch typeFlag {
	default:
		log.Fatal("")
		cmdDuplicates.Usage()
	case "fasta":
		duplicates(br)
	case "range":
		r := make(region)
		r.ReadFrom(br)
		r.WriteTo(bw)
	}
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
func duplicates(r io.Reader) (err error) {
	dupsData := make(Fasta)

	var referenceSequence, querySequence string
	var referenceLength, queryLength int64
	var referenceStart, referenceEnd, queryStart, queryEnd int64

	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		fields := strings.Fields(scanner.Text())

		if strings.HasPrefix(fields[0], ">") {
			// Header line format: >reference_contig query_contig reference_length query_length
			referenceSequence = fields[0][1:]
			querySequence = fields[1]
			if _, ok := dupsData[querySequence]; !ok {
				if queryLength, err = strconv.ParseInt(fields[3], 10, 64); err != nil {
					return err
				}
				positions := make([]byte, queryLength)
				for i := range positions {
					positions[i] = '0'
				}
				dupsData[querySequence] = positions
			}
			if _, ok := dupsData[referenceSequence]; !ok {
				if referenceLength, err = strconv.ParseInt(fields[2], 10, 64); err != nil {
					return err
				}
				positions := make([]byte, referenceLength)
				for i := range positions {
					positions[i] = '0'
				}
				dupsData[referenceSequence] = positions
			}
		} else if len(fields) == 7 {
			// Alignment line format: reference_start reference_end query_start query_end
			if referenceStart, err = strconv.ParseInt(fields[0], 10, 64); err != nil {
				return err
			}
			if referenceEnd, err = strconv.ParseInt(fields[1], 10, 64); err != nil {
				return err
			}
			if queryStart, err = strconv.ParseInt(fields[2], 10, 64); err != nil {
				return err
			}
			if queryEnd, err = strconv.ParseInt(fields[3], 10, 64); err != nil {
				return err
			}

			if isSameRegion(referenceSequence, querySequence, referenceStart, referenceEnd, queryStart, queryEnd) {
				// Skip alignments for the same region of the same contig
				break
			}

			// Swap
			if referenceStart > referenceEnd {
				referenceStart, referenceEnd = referenceEnd, referenceStart
			}
			if queryStart > queryEnd {
				queryStart, queryEnd = queryEnd, queryStart
			}

			// Mark duplicate positions
			for i := referenceStart - 1; i < referenceEnd; i++ {
				dupsData[referenceSequence][i] = '1'
			}
			for i := queryStart - 1; i < queryEnd; i++ {
				dupsData[querySequence][i] = '1'
			}
		} // else { Ignore alignment offsets }
	}

	// Print the duplicates encoded fasta
	if _, err := dupsData.WriteTo(os.Stdout); err != nil {
		return err
	}

	return nil
}

type Fasta map[string][]byte

// WriteTo implements the io.WriterTo interface.
//     >ContigName
//     <80 character wrapped contig sequence>
func (f Fasta) WriteTo(w io.Writer) (total int64, err error) {
	var contigNames []string
	for k := range f {
		contigNames = append(contigNames, k)
	}
	sort.Strings(contigNames)

	for _, contigName := range contigNames {
		positions := f[contigName]

		n, err := fmt.Fprintf(w, ">%s\n", contigName)
		total += int64(n)
		if err != nil {
			return total, err
		}

		for i := 80; i < len(positions); i += 80 {
			n, err := fmt.Fprintf(w, "%s\n", positions[i-80:i])
			total += int64(n)
			if err != nil {
				return total, err
			}
		}

		if remainder := len(positions) % 80; remainder >= 0 {
			if remainder == 0 {
				remainder = 80
			}
			n, err := fmt.Fprintf(w, "%s\n", positions[len(positions)-remainder:])
			total += int64(n)
			if err != nil {
				return total, err
			}
		}
	}

	return total, err
}

// region is a map where the key is a contig name and the value is a list of
// start and end positions for potential duplicate regions.
type region map[string][][2]int64

func (reg region) ReadFrom(r io.Reader) (nBytes int64, err error) {
	delta := mummer.Delta{}
	delta.ReadFrom(r)

	for i := range delta.Records {
		var ranges [][2]int64
		header := delta.Records[i].Header
		alignments := delta.Records[i].Alignments

		for j := range alignments {
			if !isSameRegion(header.ReferenceSequence, header.QuerySequence, alignments[j].ReferenceStart, alignments[j].ReferenceEnd, alignments[j].QueryStart, alignments[j].QueryEnd) {
				ranges = append(ranges, [2]int64{alignments[j].ReferenceStart, alignments[j].ReferenceEnd})
			}
		}

		reg[delta.Records[i].Header.ReferenceSequence] = ranges
	}

	return nBytes, nil
}

func (reg region) mergeRanges(ranges [][2]int64) [][2]int64 {
	var merge [][2]int64 = make([][2]int64, 0, len(ranges))
	var start, end int64

	for _, r := range ranges {
		if r[0] >= start && r[0] <= end && r[1] > end {
			merge[len(merge)-1][1] = r[1]
		} else if r[1] >= start && r[1] <= end && r[0] < start {
			merge[len(merge)-1][0] = r[0]
		} else {
			start, end = r[0], r[1]
			merge = append(merge, r)
		}
	}

	return merge
}

// WriteTo implements the io.WriterTo interface.
//     <contig name> <start position> <end position>
func (reg region) WriteTo(w io.Writer) (nBytes int64, err error) {
	var contigNames []string
	for key := range reg {
		contigNames = append(contigNames, key)
	}
	sort.Strings(contigNames)

	for _, contigName := range contigNames {
		for _, ranges := range reg[contigName] {
			n, err := fmt.Fprintf(w, "%s\t%d\t%d\n", contigName, ranges[0], ranges[1])
			nBytes += int64(n)
			if err != nil {
				return nBytes, err
			}
		}
	}

	return nBytes, nil
}

func isSameRegion(referenceSequence, querySequence string, referenceStart, referenceEnd, queryStart, queryEnd int64) bool {
	return referenceSequence == querySequence && (referenceStart == queryStart || referenceEnd == queryEnd)
}
