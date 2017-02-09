package matrix

import (
	"bufio"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"strconv"
	"strings"
)

type matrix struct {
	file *os.File
	buf  *bufio.Writer
}

func newMatrix(filepath string, identifiers string) *matrix {
	file, err := os.Create(filepath)
	if err != nil {
		log.Fatal(err)
	}

	buf := bufio.NewWriter(file)

	buf.Write([]byte("LocusID\tReference\t"))
	buf.WriteString(identifiers)
	buf.Write([]byte("\t#SNPcall\t#Indelcall\t#Refcall\t#CallWasMade\t#PassedDepthFilter\t#PassedProportionFilter\t#A\t#C\t#G\t#T\t#Indel\t#NXdegen\tContig\tPosition\tInDupRegion\tSampleConsensus\tCallWasMade\tPassedDepthFilter\tPassedProportionFilter\tPattern\n"))

	return &matrix{
		file: file,
		buf:  buf,
	}
}

func (m matrix) Write(line []byte) (int, error) {
	return m.buf.Write(line)
}

func (m matrix) Close() error {
	m.buf.Flush()
	return m.file.Close()
}

const (
	vcfMetadataVersions     = "##fileFormat=VCFv4.2\n##source=NASPv1.0.0\n"
	vcfMetadataDescriptions = `##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of Samples With Data">
##FILTER=<ID=NoCall,Description="No call for this sample at this position">
##FILTER=<ID=CovFail,Description="Insufficient depth of coverage for this sample at this position">
##FILTER=<ID=PropFail,Description="Insufficient proportion of reads were variant for this sample at this position">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">`
)

func newVcf(filepath string, identifiers, contigs []string) *matrix {
	file, err := os.Create(filepath)
	if err != nil {
		log.Fatal(err)
	}

	buf := bufio.NewWriter(file)

	buf.WriteString(vcfMetadataVersions)

	for i := range contigs {
		fmt.Fprintf(buf, "##contig=<ID=\"%s\",length=%d>\n", contigs[i])
	}

	for i := range identifiers {
		fmt.Fprintf(buf, "##SAMPLE=<ID=\"%[1]s\",Genomes=\"%[1]s\",Mixture=1.0>\n", identifiers[i])
	}

	buf.WriteString(vcfMetadataDescriptions)

	return &matrix{
		file: file,
		buf:  buf,
	}
}

// Pattern builds a sequence such as 11N211...
// In this example, we know:
// - The reference was called A/C/G/T because the first byte is '1'. Otherwise it would be 'N'.
// - The first sample call matches the reference because the second byte is also a '1'
// - The second sample was either: a degeneracy, not called, or failed either
//   the coverage or proportion filter because the third byte is 'N'
// - The third sample was called A/C/G/T, but was a SNP because the fourth byte is incremented.
// - The next two samples matched the reference because they are also '1'
// When a new SNP is encountered, it is assigned an incremented value.
type Pattern struct {
	next   byte
	legend map[byte]byte
}

// TODO: Would we gain anything using a sync.Pool?
func NewPattern() *Pattern {
	return &Pattern{
		next: '0',
		legend: map[byte]byte{
			'G': 0,
			'A': 0,
			'T': 0,
			'C': 0,
		},
	}
}

// Reset resets the legend for a new position sequence.
func (p *Pattern) Reset() {
	p.next = '0'
	for k := range p.legend {
		p.legend[k] = 0
	}
}

// Next returns the next byte in a position sequence.
// passFilters means the position passed the coverage and proportion filters.
func (p *Pattern) Next(call byte, passFilters bool) byte {
	if !passFilters {
		return '0'
	}

	if next, ok := p.legend[call]; !ok {
		// call != A/C/G/T
		return '0'
	} else if next != 0 {
		// `call` was previously encountered and assigned a value for this sequence.
		return next
	}

	// This is the first time `call` has been encountered.
	// Assign it an incremented value.
	p.next++
	p.legend[call] = p.next
	return p.next
}

type ContigStats struct {
	file io.Closer
	bw   *bufio.Writer
	buf  []byte
	sum  ContigStat
}

func NewContigStats(statsFolder string) (*ContigStats, error) {
	file, err := os.Create(filepath.Join(statsFolder, "general_stats.tsv"))
	if err != nil {
		return nil, err
	}
	bw := bufio.NewWriter(file)

	// Write header
	bw.Write([]byte("Contig\treference_length\treference_clean\treference_clean (%)\treference_duplicated\treference_duplicated (%)\tall_called\tall_called (%)\tall_passed_coverage\tall_passed_coverage (%)\tall_passed_proportion\tall_passed_proportion (%)\tall_passed_consensus\tall_passed_consensus (%)\tquality_breadth\tquality_breadth (%)\tany_snps\tany_snps (%)\tbest_snps\tbest_snps (%)\n\n"))

	return &ContigStats{
		file: file,
		bw:   bw,
		buf:  make([]byte, 0),
		sum: ContigStat{
			name: "Whole Genome",
		},
	}, nil
}

// Close writes the whole genome sum for all the contigs, flushes the buffer, and
// closes the file.
func (c *ContigStats) Close() error {
	c.bw.WriteByte('\n')
	c.Write(c.sum)
	c.bw.Flush()
	return c.file.Close()
}

func (c *ContigStats) Write(stat ContigStat) {
	c.sum.Add(stat)

	c.buf = c.buf[:0]

	c.buf = append(c.buf, stat.name...)
	c.buf = append(c.buf, '\t')
	c.buf = strconv.AppendInt(c.buf, int64(stat.referenceLength), 10)
	c.buf = append(c.buf, '\t')
	c.buf = strconv.AppendInt(c.buf, int64(stat.referenceGATC), 10)
	c.buf = append(c.buf, '\t')
	c.buf = strconv.AppendFloat(c.buf, float64(stat.referenceGATC)/float64(stat.referenceLength)*100, 'f', 2, 64)
	c.buf = append(c.buf, '%', '\t')
	c.buf = strconv.AppendInt(c.buf, int64(stat.referenceDuplicated), 10)
	c.buf = append(c.buf, '\t')
	c.buf = strconv.AppendFloat(c.buf, float64(stat.referenceDuplicated)/float64(stat.referenceLength)*100, 'f', 2, 64)
	c.buf = append(c.buf, '%', '\t')
	c.buf = strconv.AppendInt(c.buf, int64(stat.allCalled), 10)
	c.buf = append(c.buf, '\t')
	c.buf = strconv.AppendFloat(c.buf, float64(stat.allCalled)/float64(stat.referenceLength)*100, 'f', 2, 64)
	c.buf = append(c.buf, '%', '\t')
	c.buf = strconv.AppendInt(c.buf, int64(stat.allPassedCoverage), 10)
	c.buf = append(c.buf, '\t')
	c.buf = strconv.AppendFloat(c.buf, float64(stat.allPassedCoverage)/float64(stat.referenceLength)*100, 'f', 2, 64)
	c.buf = append(c.buf, '%', '\t')
	c.buf = strconv.AppendInt(c.buf, int64(stat.allPassedProportion), 10)
	c.buf = append(c.buf, '\t')
	c.buf = strconv.AppendFloat(c.buf, float64(stat.allPassedProportion)/float64(stat.referenceLength)*100, 'f', 2, 64)
	c.buf = append(c.buf, '%', '\t')
	c.buf = strconv.AppendInt(c.buf, int64(stat.allPassedConsensus), 10)
	c.buf = append(c.buf, '\t')
	c.buf = strconv.AppendFloat(c.buf, float64(stat.allPassedConsensus)/float64(stat.referenceLength)*100, 'f', 2, 64)
	c.buf = append(c.buf, '%', '\t')
	c.buf = strconv.AppendInt(c.buf, int64(stat.qualityBreadth), 10)
	c.buf = append(c.buf, '\t')
	c.buf = strconv.AppendFloat(c.buf, float64(stat.qualityBreadth)/float64(stat.referenceLength)*100, 'f', 2, 64)
	c.buf = append(c.buf, '%', '\t')
	c.buf = strconv.AppendInt(c.buf, int64(stat.anySnps), 10)
	c.buf = append(c.buf, '\t')
	c.buf = strconv.AppendFloat(c.buf, float64(stat.anySnps)/float64(stat.referenceLength)*100, 'f', 2, 64)
	c.buf = append(c.buf, '%', '\t')
	c.buf = strconv.AppendInt(c.buf, int64(stat.bestSnps), 10)
	c.buf = append(c.buf, '\t')
	c.buf = strconv.AppendFloat(c.buf, float64(stat.bestSnps)/float64(stat.referenceLength)*100, 'f', 2, 64)
	c.buf = append(c.buf, '%', '\n')
	c.bw.Write(c.buf)
}

type writeCloser interface {
	io.Writer
	io.Closer
}

type devNull struct{}

func (devNull) Write(p []byte) (int, error) {
	return len(p), nil
}

func (devNull) Close() error { return nil }

// Put 1 positionsPool'%', , maxBufSize positionPool
func writeMaster(chunks chan Chunk, matrixFolder, statsFolder string, names []string, isWithAllRefPos bool) {
	var withallref writeCloser
	//interruptChan := make(chan os.Signal, 1)
	//signal.Notify(interruptChan, os'%', .Interrupt)

	namesHeader := strings.Join(names, "\t")

	master := newMatrix(filepath.Join(matrixFolder, "master.tsv"), namesHeader)

	bestsnp := newMatrix(filepath.Join(matrixFolder, "bestsnp.tsv"), namesHeader)
	//bestsnpVcf := newVcf(filepath.Join(matrixFolder, "bestsnp.vcf"), names, []string{})

	missingdata := newMatrix(filepath.Join(matrixFolder, "missingdata.tsv"), namesHeader)
	//missingdataVcf := newVcf(filepath.Join(matrixFolder, "missingdata.vcf"), names, []string{})

	if isWithAllRefPos {
		withallref = writeCloser(newMatrix(filepath.Join(matrixFolder, "master_masked.tsv"), namesHeader))
	} else {
		withallref = devNull{}
	}

	// TODO: handle err
	contigStats, _ := NewContigStats(statsFolder)

	defer func() {
		log.Println("Flushing buffers, closing files, and shutting down.")
		defer master.Close()
		defer bestsnp.Close()
		//defer bestsnpVcf.Close()
		defer missingdata.Close()
		//defer missingdataVcf.Close()
		defer withallref.Close()
		defer contigStats.Close()
	}()

	outOfNSamples := []byte{'/'}
	outOfNSamples = strconv.AppendInt(outOfNSamples, int64(NUM_SAMPLES), 10)
	outOfNSamples = append(outOfNSamples, '\t')

	var prevContigName string
	var positionNum int64

	var line = line{}

	var contigStat ContigStat

	for {
		//select {
		//case <-interruptChan:
		//	fmt.Println("Caught interrupt: Flushing buffers, closing files, and shutting down.")
		//	return
		//case chunk, ok := <-chunks:
		chunk, ok := <-chunks
		if !ok {
			break
		}

		//fmt.Println("Pending", len(chunks))
		positions := <-chunk.positionsChan

		for _, position := range positions {
			// Reset the position counter when starting a new contig.
			if chunk.contigName != prevContigName {
				positionNum = 0
				prevContigName = chunk.contigName

				if contigStat.name != "" {
					contigStats.Write(contigStat)
				}

				contigStat = ContigStat{
					name: chunk.contigName,
				}
			}
			positionNum++

			contigStat.increment(position)

			// buffer the matrix line
			line.matrixFormat(outOfNSamples, chunk.contigName, positionNum, position)

			master.Write(line.buffer)

			if position.isBestSnp {
				bestsnp.Write(line.buffer)
			}

			// Replace callStringRegion with maskedCallString.
			// Only the missingdata and withallref matrices should follow.
			line.replaceCallRegion(position.maskedCallStr)

			if position.isMissingMatrix {
				missingdata.Write(line.buffer)
			}
			withallref.Write(line.buffer)

			/*
				// buffer the vcf line
				line.vcfFormat(chunk.contigName, positionNum, position)

					if position.isBestSnp {
						bestsnpVcf.Write(line.buffer)
					}

					if position.isMissingMatrix {
						missingdataVcf.Write(line.buffer)
					}
			*/

			positionPool.Put(position)
		}

		// Write final contig stat.
		positionsPool.Put(positions)
		//}
	}
	contigStats.Write(contigStat)
}

type line struct {
	startCallRegion int
	buffer          []byte
}

func (m *line) matrixFormat(outOfNSamples []byte, contigName string, positionNum int64, position *Position) []byte {
	m.buffer = m.buffer[:0]

	// LocusID
	m.buffer = append(m.buffer, contigName...)
	m.buffer = append(m.buffer, ':', ':')
	m.buffer = strconv.AppendInt(m.buffer, positionNum, 10)
	m.buffer = append(m.buffer, '\t')

	m.startCallRegion = len(m.buffer)

	for _, call := range position.callStr {
		m.buffer = append(m.buffer, call, '\t')
	}

	m.buffer = strconv.AppendInt(m.buffer, int64(position.calledSnp), 10)
	// TODO: remove static zero
	// #Indelcall
	//m.buffer = append(m.buffer, '\t', '\t')
	m.buffer = append(m.buffer, '\t', '0', '\t')
	m.buffer = strconv.AppendInt(m.buffer, int64(position.calledReference), 10)
	m.buffer = append(m.buffer, '\t')
	//
	m.buffer = strconv.AppendInt(m.buffer, int64(position.wasCalled), 10)
	m.buffer = append(m.buffer, outOfNSamples...)
	m.buffer = strconv.AppendInt(m.buffer, int64(position.passedCoverage), 10)
	m.buffer = append(m.buffer, outOfNSamples...)
	m.buffer = strconv.AppendInt(m.buffer, int64(position.passedProportion), 10)
	m.buffer = append(m.buffer, outOfNSamples...)
	//
	m.buffer = strconv.AppendInt(m.buffer, int64(position.a), 10)
	m.buffer = append(m.buffer, '\t')
	m.buffer = strconv.AppendInt(m.buffer, int64(position.c), 10)
	m.buffer = append(m.buffer, '\t')
	m.buffer = strconv.AppendInt(m.buffer, int64(position.g), 10)
	m.buffer = append(m.buffer, '\t')
	m.buffer = strconv.AppendInt(m.buffer, int64(position.t), 10)
	// TODO: remove static zero
	// #Indel
	//m.buffer = append(m.buffer, '\t', '\t')
	m.buffer = append(m.buffer, '\t', '0', '\t')
	m.buffer = strconv.AppendInt(m.buffer, int64(position.n), 10)

	m.buffer = append(m.buffer, '\t')
	m.buffer = append(m.buffer, contigName...)
	m.buffer = append(m.buffer, '\t')
	m.buffer = strconv.AppendInt(m.buffer, int64(positionNum), 10)
	m.buffer = append(m.buffer, "\t"...)
	if position.isReferenceDuplicated {
		m.buffer = append(m.buffer, "True\t"...)
	} else {
		m.buffer = append(m.buffer, "False\t"...)
	}
	if position.isAllPassConsensus {
		m.buffer = append(m.buffer, "True\t"...)
	} else {
		m.buffer = append(m.buffer, "False\t"...)
	}
	m.buffer = append(m.buffer, position.callWasMade...)
	m.buffer = append(m.buffer, '\t')
	m.buffer = append(m.buffer, position.passedDepthFilter...)
	m.buffer = append(m.buffer, '\t')
	m.buffer = append(m.buffer, position.passedProportionFilter...)
	m.buffer = append(m.buffer, '\t')
	m.buffer = append(m.buffer, position.pattern...)
	m.buffer = append(m.buffer, '\n')

	return m.buffer
}

func (l *line) vcfFormat(contigName string, positionNum int64, position *Position) []byte {
	l.buffer = l.buffer[:0]
	// #CHROM
	l.buffer = append(l.buffer, contigName...)
	l.buffer = append(l.buffer, '\t')
	// POS
	l.buffer = strconv.AppendInt(l.buffer, int64(positionNum), 10)
	// ID
	l.buffer = append(l.buffer, '\t', '.', '\t')
	// REF
	l.buffer = append(l.buffer, position.callStr[0], '\t')
	// QUAL
	l.buffer = append(l.buffer, '\t', '.', '\t')
	// FILTER
	if position.isAllCalled && position.isAllPassCoverage && position.isAllPassProportion {
		l.buffer = append(l.buffer, 'P', 'A', 'S', 'S', '\t')
	} else {
		l.buffer = append(l.buffer, 'F', 'A', 'I', 'L', '\t')
	}
	// INFO
	l.buffer = append(l.buffer, 'A', 'N', '=')
	//l.buffer = strconv.AppendInt(l.buffer, len(alts), 10)
	l.buffer = append(l.buffer, ';', 'N', 'S', '=')
	l.buffer = strconv.AppendInt(l.buffer, position.calledReference+position.calledSnp, 10)
	// FORMAT
	l.buffer = append(l.buffer, '\t', 'G', 'T', ':', 'F', 'T', '\t')
	// sample analysis columns
	//l.vcfSampleAnalysisColumns(position.pattern, position.callWasMade, position.passedDepthFilter, position.passedProportionFilter)
	return l.buffer
}

// vcfSampleAnalysisColumns appends a GT:FT vector column for each sample analysis
// and ends the lines with a \n
//     GT:FT[\tGT:FT]\n
func (l *line) vcfSampleAnalysisColumns(pattern, callWasMade, passedDepthFilter, passedProportionFilter []byte) {
	for i := range pattern {
		// Assume pattern is a string consisting of the digits 0,1,2,3,4
		if pattern[i] == '0' {
			// '0' indicates the analysis was not of sufficent quality to identify
			// Early versions of NASP used an 'N' here
			l.buffer = append(l.buffer, '.', ':')
		} else {
			// Since we know the byte is one of 1,2,3,4 we can take advantage of
			// ASCII codes to subtract 1 without parsing it as an integer.
			l.buffer = append(l.buffer, pattern[i]-1, ':')
		}

		if callWasMade[i] != 'Y' {
			l.buffer = append(l.buffer, 'N', 'o', 'C', 'a', 'l', 'l', '\t')
		} else if passedDepthFilter[i] != 'Y' {
			l.buffer = append(l.buffer, 'C', 'o', 'v', 'F', 'a', 'i', 'l', '\t')
		} else if passedProportionFilter[i] != 'Y' {
			l.buffer = append(l.buffer, 'P', 'r', 'o', 'p', 'F', 'a', 'i', 'l', '\t')
		} else {
			l.buffer = append(l.buffer, 'P', 'A', 'S', 'S', '\t')
		}
	}
	// The above loop ended the line with a \t. Since this is the last character,
	// replace it with a '\n' to terminate the line.
	l.buffer[len(l.buffer)-1] = '\n'
}

/*
def _vcf_analysis_column(pattern, analysis_stats):
    """
    Args:
        pattern (str):
        analysis_stats:

    Return:
        str: The format string for each analysis column.
    """
    for pattern_num, analysis_stat in zip(pattern, itertools.chain.from_iterable(analysis_stats)):
        try:
            gt = int(pattern_num) - 1
        except ValueError:
            # pattern_num for the analysis was 'N'
            gt = '.'
        if not analysis_stat['was_called']:
            ft = "NoCall"
        elif not analysis_stat['passed_coverage_filter']:
            ft = "CovFail"
        elif not analysis_stat['passed_proportion_filter']:
            ft = "PropFail"
        else:
            ft = "PASS"
        yield '{0}:{1}'.format(gt, ft)
*/

func (l *line) replaceCallRegion(mask []byte) {
	for i := 0; i < len(mask); i++ {
		l.buffer[2*i+l.startCallRegion] = mask[i]
	}
}
