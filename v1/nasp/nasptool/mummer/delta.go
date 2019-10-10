package mummer

import (
	"bufio"
	"errors"
	"io"
	"strconv"
	"strings"
)

var (
	// TODO: Report line number and content
	ErrUnexpectedFormat = errors.New("nasp: parse error")
)

type Header struct {
	ReferenceSequence string
	QuerySequence     string
	ReferenceLength   int64
	QueryLength       int64
}

func (h *Header) unmarshal(fields []string) (err error) {
	h.ReferenceSequence = fields[0][1:]
	h.QuerySequence = fields[1]
	if h.QueryLength, err = strconv.ParseInt(fields[3], 10, 0); err != nil {
		return err
	}
	if h.ReferenceLength, err = strconv.ParseInt(fields[2], 10, 0); err != nil {
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
type Alignment struct {
	ReferenceStart   int64
	ReferenceEnd     int64
	QueryStart       int64
	QueryEnd         int64
	NumberOfErrors   int64
	SimilarityErrors int64
	StopCodons       int64
	Distances        []int64
}

type DeltaRecord struct {
	Header     Header
	Alignments []Alignment
}

func (d *DeltaRecord) appendAlignment(fields []string) (err error) {
	var a Alignment

	if a.ReferenceStart, err = strconv.ParseInt(fields[0], 10, 0); err != nil {
		return err
	}
	if a.ReferenceEnd, err = strconv.ParseInt(fields[1], 10, 0); err != nil {
		return err
	}
	if a.QueryStart, err = strconv.ParseInt(fields[2], 10, 0); err != nil {
		return err
	}
	if a.QueryEnd, err = strconv.ParseInt(fields[3], 10, 0); err != nil {
		return err
	}
	if a.NumberOfErrors, err = strconv.ParseInt(fields[4], 10, 0); err != nil {
		return err
	}
	if a.SimilarityErrors, err = strconv.ParseInt(fields[5], 10, 0); err != nil {
		return err
	}
	if a.StopCodons, err = strconv.ParseInt(fields[6], 10, 0); err != nil {
		return err
	}

	d.Alignments = append(d.Alignments, a)

	return nil
}

func (d *DeltaRecord) appendDistance(field string) (err error) {
	distance, err := strconv.ParseInt(field, 10, 0)
	if err != nil {
		return err
	}

	lastAlignment := len(d.Alignments) - 1
	d.Alignments[lastAlignment].Distances = append(d.Alignments[lastAlignment].Distances, distance)

	return nil
}

type Delta struct {
	DataType          string
	ReferenceFilepath string
	QueryFilepath     string
	Records           []*DeltaRecord
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
// TODO: Refactor to io.ReadFrom interface
func (d *Delta) ReadFrom(r io.Reader) (n int64, err error) {
	var record *DeltaRecord

	scanner := bufio.NewScanner(r)

	// The first line is the reference and query fasta filepaths separated by a space
	if !scanner.Scan() {
		return n, ErrUnexpectedFormat
	}

	// TODO: verify fields is composed of the reference and query filepaths
	fields := strings.Fields(scanner.Text())
	d.QueryFilepath = fields[1]

	// The second line specifies the alignment data type
	if !scanner.Scan() {
		return n, ErrUnexpectedFormat
	}

	// TODO: Add PROMER?
	if d.DataType = scanner.Text(); d.DataType != "NUCMER" {
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
				d.Records = append(d.Records, record)
			}
			record = &DeltaRecord{}
			if record.Header.unmarshal(fields); err != nil {
				return n, err
			}
		case len(fields) == 7:
			// Alignment line format: reference_start reference_end query_start query_end
			if err := record.appendAlignment(fields); err != nil {
				return n, err
			}
		case len(fields) == 1:
			// Distance is a cumulative offset to the next insert/delete
			var distances []int

			distance := parseInt([]byte(fields[0]))

			if distance == 0 {
				break
			}

			for scanner.Scan() {
				distance = parseInt(scanner.Bytes())

				if distance == 0 {
					break
				} else {
					distances = append(distances, distance)
				}
			}
			if err := record.appendDistance(fields[0]); err != nil {
				return n, err
			}
		}
	}

	d.Records = append(d.Records, record)

	return n, nil
}

func parseInt(b []byte) (n int) {
	var isNegative bool

	if len(b) == 0 {
		return 0
	}

	if b[0] == '-' {
		isNegative = true
		b = b[1:]
	}

	for i := range b {
		n = n*10 + int(b[i]-'0')
	}

	if isNegative {
		return -n
	}

	return n
}
