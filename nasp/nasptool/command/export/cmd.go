package export

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strconv"

	"github.com/TGenNorth/nasp/command"
)

var cmd = &command.Command{
	UsageLine: "export [NASP matrix].tsv",
	Short:     "convert NASP matrix output to VCF format",
	Long: `
	export takes a NASP matrix output file and writes it in VCF format to a file with the same name and a .vcf extension
	`,
}

func init() {
	cmd.Run = run

	command.Register(cmd)
}

func run(cmd *command.Command, args []string) {
	/*
		defer profile.Start(&profile.Config{
			CPUProfile: true,
			//MemProfile:     true,
			//BlockProfile:   true,
			ProfilePath:    ".",
			NoShutdownHook: true,
		}).Stop()
	*/

	if len(args) == 0 {
		cmd.Usage()
	}

	for _, matrixFilepath := range args {
		if err := exportVcf(matrixFilepath); err != nil {
			log.Println(err)
			cmd.Usage()
		}
	}
}

func exportVcf(matrixFilepath string) error {
	ext := filepath.Ext(matrixFilepath)

	infile, err := os.Open(matrixFilepath)
	if err != nil {
		return err
	}
	defer infile.Close()

	outfile, err := os.Create(matrixFilepath[:len(matrixFilepath)-len(ext)] + ".vcf")
	if err != nil {
		return err
	}
	defer outfile.Close()

	br := bufio.NewReader(infile)
	bw := bufio.NewWriter(outfile)
	defer bw.Flush()

	contigs, identifiers, err := collectVcfMetadata(br)
	infile.Seek(0, os.SEEK_SET)
	br.Reset(infile)

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

// LocusID\tReference\t<sample analysis columns>\t#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\tPattern
// <contig name>::<position>
// TODO: implement io.Writer interface
func (v *vcf) Write(line []byte) (n int64, err error) {
	var start, end int

	sampleColumns, isAllPassFilters, err := v.translateSampleAnalysisColumns(line)
	if err != nil {
		return 0, err
	}

	end = bytes.IndexByte(line, ':')

	// VCF #CHROM column
	v.w.Write(line[:end])
	v.w.WriteByte('\t')

	// Starting from the end of the contig name we found earlier and adding 2
	// to skip the '::' find the position from the LocusID column
	start = end + 2
	end = start + bytes.IndexByte(line[start:], '\t')

	// Starting after the tab following the LocusID column, the callRegion is
	// the columns of nucleotide calls between the LocusID and #SNPcall columns.
	startCallRegion := end + 1
	// endCallRegion is derived starting from the Reference column counting
	// all the sample columns and delimiting tab characters.
	endCallRegion := startCallRegion + 2*v.numSamples
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
	v.w.WriteString(strconv.Itoa(len(alts)))
	v.w.WriteString(";NS=")
	start = endCallRegion + 2
	end = start + bytes.IndexByte(line[start:], '\t')
	// Copy #SNPcall column
	v.w.Write(line[start:end])

	v.w.Write(sampleColumns)

	return 0, nil
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
		return append(alts, '.')
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
		nocall   = ":NoCall\t"
		covfail  = ":CovFail\t"
		propfail = ":PropFail\t"
		pass     = ":PASS\t"
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
	// The above loop ended the line with a \t. Since this is the last character,
	// replace it with a '\n' to terminate the line.
	buf[len(buf)-1] = '\n'

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
	end := bytes.Index(buffer, []byte("\t#SNPcall"))
	if end == -1 {
		return nil, nil, errors.New("export: #SNPcall header column not found. Expected the sample identifiers to be the columns between the LocusID and #SNPcall column")
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
