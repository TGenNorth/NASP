package export

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"log"
	"os"
	"strconv"
	"strings"

	"github.com/TGenNorth/NASP/command"
)

var cmd = &command.Command{
	UsageLine: "export --type [vcf|fasta] NASP_MATRIX > output",
	Short:     "transform a NASP matrix to FASTA or VCF format",
	Long: `
transform a NASP matrix to FASTA or VCF format

The --type flag sets the output type (default: vcf)
	`,
}

var typeFlag string

func init() {
	// TODO: Update Run interface to be a function that returns an error
	cmd.Run = run
	cmd.Flag.StringVar(&typeFlag, "type", "vcf", "")
	command.Register(cmd)
}

func run(cmd *command.Command, args []string) (err error) {
	if len(args) == 0 {
		return errors.New("export: requires a NASP matrix .tsv file")
	}

	file, err := os.Open(args[0])
	if err != nil {
		return err
	}
	defer file.Close()

	switch strings.ToLower(typeFlag) {
	default:
		err = fmt.Errorf("export: unsupported output type: '%s'\n", typeFlag)
	case "vcf":
		err = exportVcf(os.Stdout, file)
	case "fasta":
		err = exportFasta(os.Stdout, file)
	}

	return err
}

func exportFasta(w io.Writer, rs io.ReadSeeker) error {
	var filePosition int64
	var lineWidth int
	referenceColumn := make([]int64, 0)
	linefeed := []byte("\n")
	br := bufio.NewReader(rs)
	bw := bufio.NewWriter(w)
	defer bw.Flush()

	line, err := br.ReadBytes('\n')
	filePosition = int64(len(line))
	if err != nil {
		return err
	}
	header := make([]byte, len(line))
	copy(header, line)
	start := bytes.Index(header, []byte("Reference")) + len("Reference") + 1
	end := bytes.Index(header, []byte("#SNPcall"))
	identifiers := bytes.Split(header[start:end-1], []byte("\t"))

	// Transform the Reference column to a fasta contig by scanning the file
	// line by line searching for the second (Reference) column.
	//
	// Track the file position of each reference call so later we can seek
	// straight to the calls for all the remaining samples.
	bw.WriteString(">Reference\n")
	for {
		line, err := br.ReadBytes('\n')
		if err != nil {
			if err != io.EOF {
				return err
			} else if len(line) == 0 {
				break
			}
		}
		// Find the reference column
		idx := bytes.IndexByte(line, '\t')
		if idx == -1 {
			// TODO: err
			log.Fatalf("%d, %s, %v\n", idx, line, err)
		}
		// The reference column is the character immediately following the first tab
		bw.Write(line[idx+1 : idx+2])
		lineWidth++
		if lineWidth == 80 {
			// Wrap the contig sequence every 80 characters
			bw.Write(linefeed)
			lineWidth = 0
		}

		referenceColumn = append(referenceColumn, filePosition+int64(idx+1))
		filePosition += int64(len(line))
	}

	// Transform each Sample Analysis column to a fasta contig.
	// The calls are found as an offset from the position of the reference calls
	p := make([]byte, 1)
	for i := range identifiers {
		lineWidth = 0
		fmt.Fprintf(bw, "\n>%s\n", identifiers[i])
		for _, filePosition := range referenceColumn {
			rs.Seek(filePosition+int64((i+1)*2), os.SEEK_SET)

			n, err := rs.Read(p)
			if n != 1 || err != nil {
				return errors.New("export: TODO err msg")
			}
			bw.Write(p)

			// Wrap sequence every 80 characters
			lineWidth++
			if lineWidth == 80 {
				if _, err := bw.Write(linefeed); err != nil {
					return err
				}
				lineWidth = 0
			}
		}
	}

	return nil
}

func exportVcf(w io.Writer, rs io.ReadSeeker) error {
	br := bufio.NewReader(rs)
	bw := bufio.NewWriter(w)
	defer bw.Flush()

	contigs, identifiers, err := collectVcfMetadata(br)
	if err != nil {
		return err
	}
	rs.Seek(0, os.SEEK_SET)
	br.Reset(rs)

	vcf := newVcf(bw, len(identifiers))
	vcf.WriteMetadataAndHeader(contigs, identifiers)

	return vcf.ReadFrom(br)
}

type vcf struct {
	w          *bufio.Writer
	numSamples int
	buffer     []byte
}

func newVcf(w *bufio.Writer, numSamples int) *vcf {
	return &vcf{
		w:          w,
		numSamples: numSamples,
		buffer:     make([]byte, 0),
	}
}

func (v *vcf) WriteMetadataAndHeader(contigs []contigmeta, identifiers [][]byte) error {
	if _, err := fmt.Fprintf(v.w, "##fileFormat=VCFv4.2\n##source=NASP\n"); err != nil {
		return err
	}

	for i := range contigs {
		if _, err := fmt.Fprintf(v.w, "##contig=<ID=\"%s\",length=%d>\n", contigs[i].name, contigs[i].length); err != nil {
			return err
		}
	}

	for i := range identifiers {
		if _, err := fmt.Fprintf(v.w, "##SAMPLE=<ID=\"%[1]s\",Genomes=\"%[1]s\",Mixture=1.0>\n", identifiers[i]); err != nil {
			return err
		}
	}

	if _, err := fmt.Fprintf(v.w, `##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##FILTER=<ID=NoCall,Description="No call for this sample at this position">
##FILTER=<ID=CovFail,Description="Insufficient depth of coverage for this sample at this position">
##FILTER=<ID=PropFail,Description="Insufficient proportion of reads were variant for this sample at this position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
`); err != nil {
		return err

	}

	if _, err := fmt.Fprintf(v.w, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", bytes.Join(identifiers, []byte("\t"))); err != nil {
		return err
	}

	return nil
}

// ALT_MONOMORPHISM is the value of the VCF ALT column when there are no SNPS
// This is expected norm, so avoid allocating for this case
var ALT_MONOMORPHISM = []byte(".")

// LocusID\tReference\t<sample analysis columns>\t#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\tPattern
// <contig name>::<position>
// TODO: implement io.Writer interface
func (v *vcf) Write(line []byte) (n int64, err error) {
	const ONE = "1"
	var start, end int
	var ns uint

	sampleColumns, isAllPassFilters, err := v.translateSampleAnalysisColumns(line)
	if err != nil {
		return 0, err
	}

	// The Matrix first column, LocusID, is the contig name and position
	// delimited by '::'. Therefore, the VCF #CHROM column is everything up
	// until the '::' delimiter
	end = bytes.IndexByte(line, ':')

	// VCF #CHROM column
	v.w.Write(line[:end])
	v.w.WriteByte('\t')

	// Starting from the end of the contig name we found earlier and skipping
	// the '::', find the position from the LocusID column
	start = end + len("::")
	end = start + bytes.IndexByte(line[start:], '\t')

	// Starting after the tab following the LocusID column, the callRegion is
	// the columns of nucleotide calls between the LocusID and #SNPcall columns.
	startCallRegion := end + len("\t")
	// endCallRegion is derived starting from the Reference column counting
	// all the sample columns and delimiting tab characters.
	endCallRegion := startCallRegion + 2*v.numSamples + 1
	alts := v.alts(line[startCallRegion:endCallRegion])

	// VCF POS column
	v.w.Write(line[start:end])

	// VCF ID column
	v.w.WriteString("\t.\t")

	// VCF REF column
	// <Reference Call>\t
	v.w.Write(line[startCallRegion : startCallRegion+2])

	v.w.Write(alts)

	// VCF QUAL FILTER FORMAT
	if isAllPassFilters {
		v.w.WriteString("\t.\tPASS\tGT:FT\t")
	} else {
		v.w.WriteString("\t.\tFAIL\tGT:FT\t")
	}

	// VCF INFO
	v.w.WriteString("AN=")
	if alts[0] == '.' {
		v.w.WriteString(ONE)
		ns = 0
	} else {
		v.w.WriteString(strconv.Itoa(len(alts) + 1))
		ns = uint(len(alts))
	}
	v.w.WriteString(";NS=")
	// Find the #SNPcall column
	start = endCallRegion + 2
	// Find the #Indelcall column
	start += bytes.IndexByte(line[start:], '\t') + 1
	// Find the #Refcall column
	start += bytes.IndexByte(line[start:], '\t') + 1
	// Find the end of #Refcall
	end = start + bytes.IndexByte(line[start:], '\t')
	ns += parseUint(line[start:end])
	v.w.WriteString(strconv.Itoa(int(ns)))

	v.w.Write(sampleColumns)

	return 0, nil
}

func parseUint(p []byte) (n uint) {
	for _, digit := range p {
		if digit >= '0' && digit <= '9' {
			n *= 10
			n += uint(digit - '0')
		} else {
			panic(fmt.Sprintf("Failed to parse '%s' as an integer"))
		}
	}
	return n
}

// alts assumes callRegion is a slice of the tab delimited call columns from a
// nasp matrix where the first call is from the Reference and the rest are from
// samples.
func (v *vcf) alts(callRegion []byte) []byte {
	alts := make([]byte, 0)
	refCall := callRegion[0]
	for i := 2; i < len(callRegion); i += 2 {
		if callRegion[i] != 'X' && callRegion[i] != 'N' && callRegion[i] != refCall && bytes.IndexByte(alts, callRegion[i]) == -1 {
			alts = append(alts, callRegion[i], ',')
		}
	}
	if len(alts) == 0 {
		return ALT_MONOMORPHISM
	}
	return alts[:len(alts)-1]
}

// translateSampleAnalysisAnalysisColumns writes
// It assumes the end of the matrix line ends with the following columns:
// - CallWasMade
// - PassedDepthFilter
// - PassedProportionFilter
// - Pattern
// It assumes the column widths are identical depending only on the number of
// samples. The Pattern column is assumed to start with one additional character
// for the reference.
// It is assumes the columns are delimited by one tab character.
func (v *vcf) translateSampleAnalysisColumns(matrixLine []byte) (buf []byte, isAllPassFilters bool, err error) {
	const (
		nocall   = ":NoCall"
		covfail  = ":CovFail"
		propfail = ":PropFail"
		pass     = ":PASS"
	)
	// Counting from the end of the line, calculate the regions of the remaining
	// matrix columns based on the number of samples.
	lineLength := len(matrixLine)
	var start, end int

	buf = v.buffer[:0]
	isAllPassFilters = true

	// -1 because the pattern column starts with a value for the reference
	start = lineLength - v.numSamples - 1
	end = lineLength
	pattern := matrixLine[start:end]

	// -1 for the tab delimiting the columns
	end = start - 1
	start = end - v.numSamples
	passedProportionFilter := matrixLine[start:end]

	// -1 for the tab delimiting the columns
	end = start - 1
	start = end - v.numSamples
	passedDepthFilter := matrixLine[start:end]

	// -1 for the tab delimiting the columns
	end = start - 1
	start = end - v.numSamples
	callWasMade := matrixLine[start:end]

	// TODO: verify all slices are the same size to avoid unexpected index out of bounds error

	for i := range callWasMade {
		// Assume pattern is a string consisting of the digits 0,1,2,3,4
		if pattern[i+1] != '0' {
			// Since we know the byte is one of 1,2,3,4 we can take advantage of
			// ASCII codes to subtract 1 without parsing it as an integer.
			buf = append(buf, '\t', pattern[i]-1)
		} else {
			// '0' indicates the analysis was not of sufficent quality to identify
			// Early versions of NASP used a pattern value of 'N'
			buf = append(buf, '\t', '.')
		}

		if callWasMade[i] == 'N' {
			buf = append(buf, nocall...)
		} else if passedDepthFilter[i] == 'N' || passedDepthFilter[i] == '?' {
			buf = append(buf, covfail...)
		} else if passedProportionFilter[i] == 'N' || passedProportionFilter[i] == '?' {
			buf = append(buf, propfail...)
		} else {
			buf = append(buf, pass...)
		}
	}
	// End the line with a linefeed
	buf = append(buf, '\n')

	return buf, isAllPassFilters, nil
}

func (v *vcf) ReadFrom(r io.Reader) error {
	scanner := bufio.NewScanner(r)

	// Skip the header line
	scanner.Scan()

	for scanner.Scan() {
		if scanner.Err() != nil {
			return scanner.Err()
		}
		if _, err := v.Write(scanner.Bytes()); err != nil {
			return err
		}
	}

	return nil
}

type contigmeta struct {
	name   string
	length int
}

// assumes reader is a NASP matrix
func collectVcfMetadata(r io.Reader) (contigs []contigmeta, identifiers [][]byte, err error) {
	var lastContigName []byte = make([]byte, 0)
	var lastContigLength int

	scanner := bufio.NewScanner(r)

	// The sample identifiers are the header columns between LocusID and #SNPcall
	scanner.Scan()
	buffer := scanner.Bytes()

	start := len("LocusID\tReference\t")
	if !bytes.Equal(buffer[:start], []byte("LocusID\tReference\t")) {
		log.Println("export: expected the first two columns of the input matrix to be LocusID and Reference; the output could be incorrect")
	}

	end := bytes.Index(buffer, []byte("\t#SNPcall"))
	if end == -1 {
		return nil, nil, errors.New("export: #SNPcall header column not found. Expected the sample identifiers to be the columns between the LocusID and #SNPcall column")
	}
	snpIndRef := []byte("\t#SNPcall\t#Indelcall\t#Refcall")
	if !bytes.Equal(buffer[end:end+len(snpIndRef)], snpIndRef) {
		log.Printf("export: expected the input matrix columns following the sample analyses to match '%s'; the output could be incorrect\n", bytes.Replace(snpIndRef, []byte("\t"), []byte(" "), -1))
	}

	// Copy the header identifiers into a new backing array so they will persist
	// after the next read.
	header := make([]byte, end-start)
	copy(header, buffer[start:end])
	identifiers = bytes.Split(header, []byte("\t"))

	for scanner.Scan() {
		buffer := scanner.Bytes()
		// Assumes each the LocusID column is <contig name>::<position>
		idx := bytes.IndexByte(buffer, ':')
		// TODO: err if idx == -1
		if !bytes.Equal(lastContigName, buffer[:idx]) {
			// Starting a new contig
			contig := contigmeta{
				name:   string(lastContigName),
				length: lastContigLength,
			}
			contigs = append(contigs, contig)
			lastContigName = lastContigName[:0]
			lastContigName = append(lastContigName, buffer[:idx]...)
			lastContigLength = 0
		}
		lastContigLength++
	}

	contig := contigmeta{
		name:   string(lastContigName),
		length: lastContigLength,
	}

	contigs = append(contigs, contig)

	return contigs, identifiers, err
}
