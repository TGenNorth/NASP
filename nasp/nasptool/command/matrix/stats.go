package matrix

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"path"
	"strconv"
	"strings"
)

type Contigs struct {
	Reference  []byte
	Duplicates []byte
	Analyses   [][]byte
}

type ContigStat struct {
	name                string
	referenceLength     int
	referenceGATC       int
	referenceDuplicated int
	allCalled           int
	allPassedCoverage   int
	allPassedProportion int
	allPassedConsensus  int
	qualityBreadth      int
	anySnps             int
	bestSnps            int
}

func NewContigStat(name string, minCoverage int, minProportion float64) *ContigStat {
	return &ContigStat{
		name: name,
	}
}

func (c *ContigStat) increment(position *Position) {
	c.referenceLength++

	if position.isReferenceGATC {
		c.referenceGATC++
	}
	if position.isReferenceDuplicated {
		c.referenceDuplicated++
	}
	if position.isAllCalled {
		c.allCalled++
	}
	if position.isAllPassCoverage {
		c.allPassedCoverage++
	}
	if position.isAllPassProportion {
		c.allPassedProportion++
	}
	if position.isAllPassConsensus {
		c.allPassedConsensus++
	}
	if position.isAllQualityBreadth {
		c.qualityBreadth++
	}
	if position.isAnySnp {
		c.anySnps++
	}
	if position.isBestSnp {
		c.bestSnps++
	}
}

func (c *ContigStat) Add(other ContigStat) {
	c.referenceLength += other.referenceLength
	c.referenceGATC += other.referenceGATC
	c.referenceDuplicated += other.referenceDuplicated
	c.allCalled += other.allCalled
	c.allPassedCoverage += other.allPassedCoverage
	c.allPassedProportion += other.allPassedProportion
	c.allPassedConsensus += other.allPassedConsensus
	c.qualityBreadth += other.qualityBreadth
	c.anySnps += other.anySnps
	c.bestSnps += other.bestSnps
}

type SampleStat struct {
	//name                   string
	//aligner                string
	//snpcaller              string
	wasCalled              int64
	passedCoverageFilter   int64
	passedProportionFilter int64
	qualityBreadth         int64
	calledReference        int64
	calledSnp              int64
	calledDegen            int64
}

func (s *SampleStat) Increment(call, coverage, proportion, quality, ref, snp, degen bool) {
	if call {
		s.wasCalled++
	}
	if coverage {
		s.passedCoverageFilter++
	}
	if proportion {
		s.passedProportionFilter++
	}
	if quality {
		s.qualityBreadth++
	}
	if ref {
		s.calledReference++
	}
	if snp {
		s.calledSnp++
	}
	if degen {
		s.calledDegen++
	}
}

type SampleStats []SampleStat

var refLength int64

func (s SampleStats) Aggregate(ch chan SampleStats) {
	for stats := range ch {
		for i := range s {
			s[i].wasCalled += stats[i].wasCalled
			s[i].passedCoverageFilter += stats[i].passedCoverageFilter
			s[i].passedProportionFilter += stats[i].passedProportionFilter
			s[i].qualityBreadth += stats[i].qualityBreadth
			s[i].calledReference += stats[i].calledReference
			s[i].calledSnp += stats[i].calledSnp
			s[i].calledDegen += stats[i].calledDegen
		}
		statsPool.Put(stats)
	}
}

func (s SampleStats) Reset() {
	for i := range s {
		s[i].wasCalled = 0
		s[i].passedCoverageFilter = 0
		s[i].passedProportionFilter = 0
		s[i].qualityBreadth = 0
		s[i].calledReference = 0
		s[i].calledSnp = 0
		s[i].calledDegen = 0
	}
}

func (s SampleStats) WriteStats(identifiers []string, statsFolder string) {
	var buf []byte

	fmt.Printf("WriteStats: %s\n", path.Join(statsFolder, "sample_stats.tsv"))
	file, err := os.Create(path.Join(statsFolder, "sample_stats.tsv"))
	if err != nil {
		log.Fatal(err)
	}
	bw := bufio.NewWriter(file)
	defer bw.Flush()
	bw.Write([]byte("Sample\tSample::Analysis\twas_called\twas_called (%)\tfailed_coverage_filter\tfailed_coverage_filter (%)\tfailed_proportion_filter\tfailed_proportion_filter (%)\tquality_breadth\tquality_breadth (%)\tcalled_reference\tcalled_reference (%)\tcalled_snp\tcalled_snp (%)\tcalled_degen\tcalled_degen (%)\n\n"))

	// Whole Genome any summary
	buf = s.statLine(buf, s[len(s)-2])
	bw.WriteString("[any]\t\t")
	bw.Write(buf)

	// Whole Genome all summary
	buf = s.statLine(buf, s[len(s)-1])
	bw.WriteString("[all]\t\t")
	bw.Write(buf)
	bw.WriteByte('\n')

	// Stats for each
	for i := 0; i < len(s)-2; i++ {
		buf = s.statLine(buf, s[i])

		idx := strings.IndexByte(identifiers[i], ':')
		if idx == -1 {
			// TODO:
			//log.Printf("DEBUG: expected identifier '%s' to contain a colon\n", identifiers[i])
			idx = len(identifiers[i])
		}
		sampleName := identifiers[i][:idx]
		bw.WriteString(sampleName)
		bw.WriteByte('\t')
		bw.WriteString("[any]\t")
		bw.Write(buf)

		bw.WriteString(sampleName)
		bw.WriteByte('\t')
		bw.WriteString("[all]\t")
		bw.Write(buf)

		bw.WriteString(sampleName)
		bw.WriteByte('\t')
		bw.WriteString(identifiers[i])
		bw.WriteByte('\t')
		bw.Write(buf)

		bw.WriteByte('\n')
	}
}

func (s SampleStats) statLine(buf []byte, stat SampleStat) []byte {
	// Jason Sahl, in the 2015-11-03 informatics meeting, suggested the sample
	// stats coverage/proportion filters should show the number of positions
	// that failed the filter instead of the number that passed. In addition,
	// the number should be relative to the number of called positions.
	//
	// Because fastas contain no coverage/proportion information, they are
	// treated as if all positions pass coverage/proportion. Since all positions
	// pass coverage/proportion, but not all positions are called the stats file
	// can show a fasta with a negative number of failed coverage/proportion
	// positions.
	//
	// The counters are set to zero to handle this scenario.
	var failedCoveragePercent, failedProportionPercent float64
	failedCoverage := stat.wasCalled - stat.passedCoverageFilter
	failedProportion := stat.wasCalled - stat.passedProportionFilter
	if failedCoverage < 0 {
		failedCoverage = 0
	}
	if failedProportion < 0 {
		failedProportion = 0
	}

	if stat.wasCalled == 0 {
		// NAN guard
		failedCoveragePercent = 0.0
		failedProportionPercent = 0.0
	} else {
		failedCoveragePercent = float64(failedCoverage) / float64(stat.wasCalled)
		failedProportionPercent = float64(failedCoverage) / float64(stat.wasCalled)
	}

	buf = buf[:0]
	buf = strconv.AppendInt(buf, stat.wasCalled, 10)
	buf = append(buf, '\t')
	buf = strconv.AppendFloat(buf, float64(stat.wasCalled)/float64(refLength)*100, 'f', 2, 64)
	buf = append(buf, '%', '\t')
	buf = strconv.AppendInt(buf, failedCoverage, 10)
	buf = append(buf, '\t')
	buf = strconv.AppendFloat(buf, failedCoveragePercent, 'f', 2, 64)
	buf = append(buf, '%', '\t')
	buf = strconv.AppendInt(buf, failedProportion, 10)
	buf = append(buf, '\t')
	buf = strconv.AppendFloat(buf, failedProportionPercent, 'f', 2, 64)
	buf = append(buf, '%', '\t', '\t') // Blank column
	buf = strconv.AppendInt(buf, refLength, 10)
	buf = append(buf, '\t')
	buf = strconv.AppendInt(buf, stat.qualityBreadth, 10)
	buf = append(buf, '\t')
	buf = strconv.AppendFloat(buf, float64(stat.qualityBreadth)/float64(refLength)*100, 'f', 2, 64)
	buf = append(buf, '%', '\t')
	buf = strconv.AppendInt(buf, stat.calledReference, 10)
	buf = append(buf, '\t')
	buf = strconv.AppendFloat(buf, float64(stat.calledReference)/float64(refLength)*100, 'f', 2, 64)
	buf = append(buf, '%', '\t')
	buf = strconv.AppendInt(buf, stat.calledSnp, 10)
	buf = append(buf, '\t')
	buf = strconv.AppendFloat(buf, float64(stat.calledSnp)/float64(refLength)*100, 'f', 2, 64)
	buf = append(buf, '%', '\t')
	buf = strconv.AppendInt(buf, stat.calledDegen, 10)
	buf = append(buf, '\t')
	buf = strconv.AppendFloat(buf, float64(stat.calledDegen)/float64(refLength)*100, 'f', 2, 64)
	buf = append(buf, '%', '\n')

	return buf
}
