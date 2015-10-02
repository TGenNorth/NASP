package matrix

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"strconv"
)

var (
	// Used to identify the Header row
	chromColumnHeader = []byte("#CHROM")
	emptyVcfCalls     = bytes.Repeat([]byte("X"), defaultBufSize)
	emptyVcfCovs      = make([]int, defaultBufSize)
	emptyVcfProps     = make([]float64, defaultBufSize)
	ar                = []byte(";AR=")
	dp4               = []byte(";DP4=")
)

type aligner int

const (
	gatk aligner = iota
	varscan
	samtools
	solsnp
	unknown
)

type Vcf struct {
	aligner    aligner
	identifier string
	index      map[string][2]int64
	// TODO: Replace with ReadSeeker interface
	rd          *os.File
	br          *bufio.Reader
	calls       []byte
	coverages   []int
	proportions []float64
	position    int
	remaining   int64
	// The last record read if its position beyond the buffer.
	record vcfRecord
}

func NewVcf(naspFile NaspFile) (*Vcf, error) {
	file, err := os.Open(naspFile.Filepath)
	if err != nil {
		return nil, fmt.Errorf("NewVcf(\"%#v\") %s", naspFile, err.Error())
	}

	aligner := unknown
	switch naspFile.Aligner {
	case "GATK":
		aligner = gatk
	case "Samtools":
		aligner = samtools
	case "VarScan":
		aligner = varscan
	case "SolSnp":
		aligner = solsnp
	}

	v := &Vcf{
		aligner:    aligner,
		identifier: fmt.Sprintf("%s::%s,%s", naspFile.Name, naspFile.Aligner, naspFile.Snpcaller),
		index:      make(map[string][2]int64),
		rd:         file,
		br:         bufio.NewReader(file),
		//buf:         bytes.NewBuffer(make([]byte, 0, defaultBufSize)),
		calls:       make([]byte, defaultBufSize),
		coverages:   make([]int, defaultBufSize),
		proportions: make([]float64, defaultBufSize),
		record:      make(vcfRecord, 10),
	}
	v.indexContigs()

	return v, nil
}

// TODO: refactor to Identifier
func (v *Vcf) Name() string {
	return v.identifier
}

func (v *Vcf) ReadPositions(n int) ([]byte, []int, []float64, error) {
	//fmt.Printf("Vcf.ReadPositions(%d)\n", n)
	var err error
	var line []byte
	var offset int

	if v.position == 0 {
		return emptyVcfCalls[:n], emptyVcfCovs[:n], emptyVcfProps[:n], nil
	}

	recordPosition := v.record.Position()

	for currentPosition, lastPosition := v.position, v.position+n; currentPosition < lastPosition && lastPosition >= recordPosition && v.remaining >= 0; {
		offset = currentPosition - v.position

		if currentPosition == recordPosition {
			call := v.record.Call()
			coverage := v.record.Coverage()
			v.calls[offset] = call
			v.coverages[offset] = coverage
			v.proportions[offset] = v.record.Proportion(call, coverage, v.aligner)

			currentPosition++
			offset++
		} else if currentPosition > recordPosition {
			// Skip duplicate positions
			/*
				if *keepLastDuplicate {
					if currentPosition == v.position {
						log.Printf("FIXME: Duplicate was on the border between chunks. Unable to change prev value for %s position %d\n", v.identifier, recordPosition)
					} else if currentPosition == recordPosition+1 {
						call := v.record.Call()
						coverage := v.record.Coverage()
						v.calls[offset-1] = call
						v.coverages[offset-1] = coverage
						v.proportions[offset-1] = v.record.Proportion(call, coverage, v.aligner)
					} else {
						log.Fatal("TODO: This is a debugging error. Expected the record position to only be behind the current position if the position was listed more than once (GATK indels)")
					}
				}
			*/
		} else {
			gap := recordPosition - currentPosition
			// The next record is within the buffer, but there are empty
			// positions between here and there.
			copy(v.calls[offset:offset+gap], emptyVcfCalls)
			copy(v.coverages[offset:offset+gap], emptyVcfCovs)
			copy(v.proportions[offset:offset+gap], emptyVcfProps)

			currentPosition = recordPosition
			continue
			// FIXME: if recordPosition fails to parse we could have an infinite loop
		}

		// Read the next record
		if line, err = v.br.ReadSlice('\n'); err != nil {
			break
		}
		v.record.Split(line)
		v.remaining--
		recordPosition = v.record.Position()
	}

	if offset < n {
		// Fill remaining positions with empty values.
		copy(v.calls[offset:n], emptyVcfCalls)
		copy(v.coverages[offset:n], emptyVcfCovs)
		copy(v.proportions[offset:n], emptyVcfProps)
	}

	v.position += n

	return v.calls[:n], v.coverages[:n], v.proportions[:n], nil
}

func (v *Vcf) SeekContig(name string) (err error) {
	if index, ok := v.index[name]; ok {
		v.position = 1
		v.remaining = index[1] - 1
		// Seek to the filePosition portion of the contig index.
		_, err = v.rd.Seek(index[0], os.SEEK_SET)
		v.br.Reset(v.rd)
		//v.buf.Reset()
		// TODO: Launch goroutine to parse records / yield position chunks
		line, err := v.br.ReadSlice('\n')
		if err != nil {
			return err
		}
		v.record.Split(line)
	} else {
		// 0 is a special value indicating the contig is not present.
		v.position = 0
	}
	return err
}

func (v *Vcf) indexContigs() error {
	var err error
	// TODO: Rename priorContigLength -> priorContigRows
	var priorContigLength, priorFilePosition, filePosition int64
	var line, contigName []byte
	priorContigName := make([]byte, 0, 10)

	// Skip metadata
	for {
		line, err = v.br.ReadSlice('\n')
		filePosition += int64(len(line))
		if err == nil {
			if line, err = v.br.Peek(6); err != nil {
				// TODO: malformed vcf no header
				panic(err)
			}
			if bytes.Equal(line, chromColumnHeader) {
				// Found the header line
				break
			}
		} else if err != bufio.ErrBufferFull {
			// NOTE: bufio.ErrBufferFull is generally ok. Complete lines are not
			// required until after the metadata. The header row won't be
			// detected if the error occurs on that line.
			// TODO: malformed vcf no header
			panic(fmt.Errorf("%s: %s\n", v.identifier, err.Error()))
		}
	}

	// Skip header
	// TODO: Get sample name and Error Handling
	line, err = v.br.ReadSlice('\n')
	if err != nil {
		return err
	}
	filePosition += int64(len(line))

	// NOTE: bufio.ErrBufferFull not handled. Could fail if the line is too long.

	// When the indexing loop encounters a new contig or EOF, it stores the
	// name, filePosition, and length of the prior contig in the index.
	// The first contig is a special case, because there is no prior contig.
	//
	// At this point, `line` contains the mandatory header row and filePosition
	// is on the first contig.
	priorFilePosition = filePosition
	if line, err = v.br.ReadSlice('\t'); err != nil {
		// TODO: The file should not be empty.
		// The file should have at least one contig position.
		panic(fmt.Errorf("%s: %s\n", v.identifier, err.Error()))
	}
	priorContigName = append(priorContigName, line...)
	filePosition += int64(len(priorContigName))
	if line, err = v.br.ReadSlice('\n'); err != nil {
		// TODO: The file should not be empty
		// The contig position should contain more than the contig name column.
		panic(fmt.Errorf("%s: %s\n", v.identifier, err.Error()))
	}
	filePosition += int64(len(line))
	priorContigLength = 1

	//tick := time.Tick(1 * time.Second)
	//var fileSize int64 = -1
	//if fi, err := v.rd.Stat(); err != nil {
	//	fileSize = fi.Size()
	//}

	for {
		//select {
		//case <-tick:
		//	fmt.Printf("\t%s: %d / %d indexed\n", v.name, filePosition, fileSize)
		//default:
		//}
		if contigName, err = v.br.ReadSlice('\t'); err != nil && err != io.EOF {
			return err
		}

		// Found the file position of a new contig
		if !bytes.Equal(contigName, priorContigName) || err == io.EOF {
			// Discard the trailing \t and cast to string to create an
			// immutable copy that won't be lost on the next call to ReadSlice
			name := string(priorContigName[:len(priorContigName)-1])
			v.index[name] = [2]int64{priorFilePosition, priorContigLength}
			priorContigName = priorContigName[:0]
			priorContigName = append(priorContigName, contigName...)
			//priorContigName = contigName
			priorFilePosition = filePosition
			priorContigLength = 0
		}

		// Wait until after we know this is not a new contig before advancing
		// the filePosition or priorContigLength counter.
		filePosition += int64(len(contigName))
		priorContigLength++

		if line, err = v.br.ReadSlice('\n'); err != nil {
			return err
		}
		filePosition += int64(len(line))
	}
}

const (
	chrom = iota
	pos
	id
	ref
	alt
	qual
	filter
	info
	format
	sample
)

type vcfRecord [][]byte

func NewVcfRecord() vcfRecord {
	return make([][]byte, 10)
}

// Split slices a VCF line by column.
// It assumes: #CHROM POS ID REF ALT QUAL FILTER INFO FORMAT <SAMPLE>
func (v vcfRecord) Split(s []byte) {
	start := 0
	na := 0
	for i := 0; i+1 <= len(s) && na+1 < 10; i++ {
		if s[i] == '\t' {
			v[na] = s[start:i]
			na++
			start = i + 1
		}
	}
	v[na] = s[start:]
}

func (v vcfRecord) Call() (call byte) {
	if len(v[alt]) == 1 && v[alt][0] == '.' {
		// Sample call matches reference
		return v[ref][0]
	}

	idx := v.gtIndex()
	switch idx {
	default:
		alts := split(v[alt], ',')
		// Handle a scenario such as this from a Varscan output:
		// +-----+-----+
		// | REF | ALT |
		// +-----+-----+
		// | GTT | TT  |
		// +-----+-----+
		if len(v[ref]) > 1 && bytes.Equal(v[ref][1:], alts[idx-1]) {
			call = v[ref][0]
		} else {
			call = alts[idx-1][0]
		}
	case -1:
		// GT is undefined or '.'
		call = 'X'
	case 0:
		call = v[ref][0]
	}
	return call
}

var (
	infoDp   = []byte("DP=")
	formatDp = []byte(":DP")
)

func (v vcfRecord) Coverage() int {
	info := v[info]
	coverage := -1

	// Search the SAMPLE column for the DP sub-field
	if idx := v.findSampleSubfield(formatDp); idx != -1 {
		return parseInt(v[sample][idx:])
	}

	// Search the INFO column for the DP or ADP sub-field.
	for i := 0; i < len(info); {
		offset := bytes.Index(info[i:], infoDp)

		if offset < 0 {
			// Field not found.
			break
		}

		i += offset

		// Verify this really is the DP or ADP sub-field, not one that
		// happens to end with 'DP'.
		if offset == 0 || info[i-1] == ';' || info[i-1] == 'A' && (offset == 1 || info[i-2] == ';') {
			// offset == 0: Found at the start of column. '\t' was discarded when the line was split into columns.
			coverage = parseInt(info[i+len(infoDp):])
			break
		} else {
			// Field happens to end with DP, move along.
			i += len(infoDp)
		}
	}

	return coverage
}

func (v vcfRecord) Proportion(call byte, coverage int, t aligner) float64 {
	isSnp := v[ref][0] != call

	switch t {
	case gatk, varscan:
		if proportion := v.gatkVarscan(isSnp, coverage); proportion != -1.0 {
			return proportion
		}
	case samtools:
		if proportion := v.samtools(isSnp, coverage); proportion != -1.0 {
			return proportion
		}
	case solsnp:
		if proportion := v.solSnp(isSnp); proportion != -1.0 {
			return proportion
		}
	default:
		if proportion := v.gatkVarscan(isSnp, coverage); proportion != -1.0 {
			return proportion
		} else if proportion := v.solSnp(isSnp); proportion != -1.0 {
			return proportion
		} else if proportion := v.samtools(isSnp, coverage); proportion != -1.0 {
			return proportion
		}
	}

	// if proportion := v.gatkVarscan(isSnp, coverage); proportion != -1.0 {
	// 	return proportion
	// } else if proportion := v.solSnp(isSnp); proportion != -1.0 {
	// 	return proportion
	// } else if proportion := v.samtools(isSnp, coverage); proportion != -1.0 {
	// 	return proportion
	// }
	return -1.0
}

var ad = []byte(":AD")
var rd = []byte(":RD")

func (v vcfRecord) gatkVarscan(isSnp bool, coverage int) float64 {
	if coverage < 1 || len(v[format]) < 1 {
		return -1.0
	}

	idx := v.findSampleSubfield(ad)
	if idx == -1 {
		return -1.0
	}

	callDepths := v[sample][idx:]
	if endCallDepths := bytes.IndexByte(callDepths, ':'); endCallDepths != -1 {
		callDepths = callDepths[:endCallDepths]
	}

	// GATK
	//if start, gt := 0, v.gtIndex(); gt > -1 {
	start, gt := 0, v.gtIndex()
	if gt > -1 && count(callDepths, ',') > 0 {
		for i := 0; i < gt; i++ {
			idx := bytes.IndexByte(callDepths, ',')
			if idx == -1 {
				//
				start = -1
				break
			}
			start += idx + 1
		}
		if start != -1 {
			return float64(parseInt(callDepths[start:])) / float64(coverage)
		}
	}

	// VarScan
	if isSnp {
		return float64(parseInt(callDepths)) / float64(coverage)
	}

	if idx := v.findSampleSubfield(rd); idx != -1 {
		return float64(parseInt(v[sample][idx:])) / float64(coverage)
	}

	return -1
}

// FIXME?: findSampleSubfield assumes subfield is a string such as ":AD"
// Remove requirement the subfield starts with a ':'
func (v vcfRecord) findSampleSubfield(subfield []byte) int {
	var nField int

	// If FORMAT is empty, the loop will never run and coverage will be 0
	// instead of -1 indicating not found.
	if len(v[format]) < 1 {
		return -1
	}

	// Find the subfield in the FORMAT column.
	for start, i := 0, 0; start < len(v[format]); start += i + len(subfield) {
		sub := v[format][start:]
		i = bytes.Index(sub, subfield)
		if i == -1 {
			if bytes.HasPrefix(v[format], subfield[1:]) && (len(v[format]) == len(subfield[1:]) || v[format][len(subfield)] == ':') {
				// This should be an unusual case. The subfield is the first
				// field and does not happen to be the prefix of another field.
				break
			} else {
				// The subfield is undefined
				return -1
			}
		}

		// Make sure we matched the correct field and not one that happens to
		// end with the same characters.
		// Either is it the last subfield or is followed by another.
		if len(sub[i:]) == len(subfield) || sub[i+len(subfield)] == ':' {
			// +1 to be sure the subfield ':' is included in the count.
			// +1 could just as easily be len(subfield).
			nField = count(v[format][:start+i+1], ':')
			break
		}
	}

	var start int

	// Find the corresponding subfield in the SAMPLE column.
	for i := 0; i < nField; i++ {
		idx := bytes.IndexByte(v[sample][start:], ':')
		if idx == -1 {
			// TODO: replace with error?
			panic(fmt.Errorf("Expected to find subfield %d in the SAMPLE column corresponding to the FORMAT column %s subfield.\nRecord: %s", nField, subfield, v))
		}
		start += idx + 1
	}

	return start
}

// samtools returns a value based on the INFO DP4 column or -1.
func (v vcfRecord) samtools(isSnp bool, coverage int) float64 {
	if coverage < 1 {
		return -1
	}

	info := v[info]
	start := 0

	if len(info) < len(dp4) {
		return -1.0
	}

	if idx := bytes.Index(info, dp4); idx != -1 {
		// DP4 is somewhere in the INFO column
		start = idx + len(dp4)
	} else if bytes.Equal(info[:len(dp4)-1], dp4[1:]) {
		// DP4 is the first sub-field.
		start = len(dp4) - 1
	} else {
		// DP4 is undefined.
		return -1
	}

	var n1, n2 int
	// DP4 is a comma separated list of numbers.
	if isSnp {
		// Parse the last 2 numbers as n1 and n2.
		var idx int

		if idx = bytes.IndexByte(info[start:], ','); idx == -1 {
			// TODO: log warning
			return -1
		}
		start += idx + 1
		if idx = bytes.IndexByte(info[start:], ','); idx == -1 {
			// TODO: log warning
			return -1
		}
		start += idx + 1

		n1 = parseInt(info[start:])

		if idx = bytes.IndexByte(info[start:], ','); idx == -1 {
			// TODO: log warning
			return -1
		}
		n2 = parseInt(info[start+idx+1:])

		return float64(n1+n2) / float64(coverage)
	} else {
		// Parse the first 2 numbers as n1 and n2.
		if n1 = parseInt(info[start:]); n1 == -1 {
			// TODO: log warning
			return -1
		}
		if idx := bytes.IndexByte(info, ','); idx == -1 {
			// TODO: log warning
			return -1
		} else if n2 = parseInt(info[idx+1:]); n2 == -1 {
			// TODO: log warning
			return -1
		}

		return float64(n1+n2) / float64(coverage)
	}
}

// solSnp returns the value of the INFO column AR sub-field else -1 if undefined
// or an error occurs while parsing.
func (v vcfRecord) solSnp(isSnp bool) float64 {
	info := v[info]
	// start and end is the index range for the digits of the AR sub-field.
	start, end := 0, 0

	if len(info) < len(ar) {
		return -1.0
	}

	if bytes.Equal(info[:len(ar)-1], ar[1:]) {
		// AR is the first sub-field.
		start = len(ar) - 1
	} else if idx := bytes.Index(info, ar); idx != -1 {
		// AR is preceeded by at least one sub-field.
		start = idx + len(ar)
	} else {
		return -1
	}

	if idx := bytes.IndexByte(info[start:], ';'); idx != -1 {
		// AR is followed by at least one sub-field.
		end = idx
	} else {
		// AR is the last sub-field.
		end = len(info[start:])
	}

	if f, err := strconv.ParseFloat(string(info[start:end]), 64); err != nil {
		// TODO: log warning
		return -1
	} else {
		return f
	}
}

// gtIndex returns the first index value from the GT field or -1 if '.' or undefined.
func (v vcfRecord) gtIndex() (n int) {
	// VCFv4.2 Section 1.4.2 Genotype fields:
	// "The first sub-field must always be the genotype (GT) if it is present."
	if len(v[format]) > 1 && v[format][0] == 'G' && v[format][1] == 'T' {
		b := v[sample]

		if b[0] == '.' {
			return -1
		}

		n = parseInt(b)
	}

	return n
}

func (v vcfRecord) Position() int {
	return parseInt(v[pos])
}

func parseInt(b []byte) (n int) {
	for _, v := range b {
		// Stop when the first non-numeric character is encounter.
		// This includes the sub-field delimiter ':',
		// undefined value '.' as well as '|' or '/'
		// indicating phased or unphased diploid and haploid calls.
		if '0' > v || v > '9' {
			break
		}
		n *= 10
		n += int(v - '0')
	}

	return n
}

// split is a simplified version of the stdlib bytes.Split()
func split(s []byte, sep byte) [][]byte {
	n := count(s, sep) + 1
	start := 0
	a := make([][]byte, n)
	na := 0
	for i := 0; i+1 <= len(s) && na+1 < n; i++ {
		if s[i] == sep {
			a[na] = s[start:i]
			na++
			start = i + 1
		}
	}
	a[na] = s[start:]
	return a[0 : na+1]
}

// count is a simplified version of the stdlib bytes.Count()
func count(s []byte, sep byte) (n int) {
	for i := 0; i < len(s); n++ {
		if o := bytes.IndexByte(s[i:], sep); o < 0 {
			break
		} else {
			i += o + 1
		}
	}

	return n
}
