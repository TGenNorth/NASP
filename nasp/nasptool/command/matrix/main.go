package matrix

import (
	"fmt"
	"io"
	"log"
	"runtime"
	"sort"
	"strings"
	"sync"
	"time"

	//"github.com/davecheney/profile"
)

const (
	// defaultBufSize is the maximum number of calls read from a contig at a time.
	// defaultBufSize * numFiles * maxQueuedBlocks roughly translates to RAM usage.
	defaultBufSize  = 4096
	maxQueuedBlocks = 20
)

var NUM_SAMPLES int
var t0 = time.Now()

/*
var refPath = flag.String("reference-fasta", "", "Path to the reference.fasta against which samples are compared")
var dupPath = flag.String("reference-dups", "", "Path to the duplicates.txt file marking duplicated positions")
var matrixPath = flag.String("matrix-folder", "", "Path to the output folder for matrices")
var statsPath = flag.String("stats-folder", "", "Path to the output folder for statistics")
var minCoverage = flag.Int("minimum-coverage", -1, "Filter positions below this coverage/depth threshold")
var minProportion = flag.Float64("minimum-proportion", -1.0, "Filter positions below this proportion threshold")
var dtoFile = flag.String("dto-file", "", "Path to the matrix_dto.xml file")
var numThreads = flag.Int("num-threads", runtime.NumCPU(), "Max number of CPUs that can be executing simultaneously")
var _ = flag.String("mode", "", "NoOp - deprecated flag for backwards compatibility")

var keepLastDuplicate = flag.Bool("keep-last", false, "If a VCF has multiple records for the same position, keep the last one")
var manualSort = flag.Bool("manual-sort", false, "Debug for developer")

func main() {
	flag.Parse()

	Run(*numThreads, *dtoFile, *refPath, *dupPath, *statsPath, *matrixPath, *minCoverage, *minProportion, flag.Args()...)
}
*/

func Run(numThreads int, dtoFile, refPath, dupPath, statsPath, matrixPath string, minCoverage int, minProportion float64, files ...string) error {
	var wg sync.WaitGroup
	runtime.GOMAXPROCS(numThreads)

	// Begin Development Profiling
	/*
		defer profile.Start(&profile.Config{
			CPUProfile: true,
			//MemProfile:     true,
			//BlockProfile:   true,
			ProfilePath:    ".",
			NoShutdownHook: true,
		}).Stop()
	*/

	t0 := time.Now()
	defer func() {
		log.Println(time.Now().Sub(t0))
	}()
	// End Development Profiling

	dto, err := NewDto(dtoFile, refPath, dupPath, matrixPath, statsPath, minCoverage, minProportion, files)
	if err != nil {
		return err
	}

	analyses := NewSampleAnalyses(dto.AllFiles...)
	sort.Sort(analyses)
	log.Println("Samples indexed", time.Now().Sub(t0))

	names := make([]string, len(analyses))
	for i := range analyses {
		names[i] = analyses[i].Name()
	}

	NUM_SAMPLES = len(analyses)

	reference, err := NewReference(dto.Reference, dto.Duplicates)
	if err != nil {
		return err
	}

	// Queue limits the rate positions are read from the files and ensures the
	// chunks are assembled in the correct order.
	queue := make(chan Chunk, maxQueuedBlocks)
	statsChan := make(chan SampleStats, maxQueuedBlocks)
	defer func() {
		wg.Add(1)
		close(queue)
		wg.Wait()
		wg.Add(1)
		close(statsChan)
		wg.Wait()
	}()

	go func(queue chan Chunk, matrixFolder, statsFolder string, withAllRefPos bool, names []string) {
		defer wg.Done()
		// queue must be closed to exit
		writeMaster(queue, matrixFolder, statsFolder, names, withAllRefPos)
	}(queue, dto.MatrixFolder, dto.StatsFolder, dto.FilterMatrixFormat == "include_allref_pos", names)

	go func(statsChan chan SampleStats, names []string, statsFolder string) {
		defer wg.Done()
		// statsCh must be closed to exit
		stats := statsPool.Get().(SampleStats)
		stats.Aggregate(statsChan)
		stats.WriteStats(names, statsFolder)
	}(statsChan, names, dto.StatsFolder)

	return readPositions(queue, statsChan, reference, analyses, dto.MinCoverage, dto.MinProportion)
}

type Chunk struct {
	contigName    string
	positionsChan chan []*Position
}

// Read positions from the reference and sample files
func readPositions(queue chan Chunk, statsChan chan SampleStats, reference *Reference, analyses SampleAnalyses, minCoverage int, minProportion float64) (err error) {
	var name string
	var ref, dup []byte
	var sampleChunk *SampleChunk

	for {
		if name, err = reference.NextContig(); err == io.EOF {
			err = nil
			break
		} else if err != nil {
			return err
		}

		if err = analyses.SeekContig(name); err != nil {
			return err
		}

		log.Printf("Scanning %s %s\n", name, time.Now().Sub(t0))

		for isPrefix := true; isPrefix; {
			ref, dup, isPrefix, err = reference.ReadPositions(defaultBufSize)
			if err != nil {
				// TODO return a type that implements the error interface instead of a raw error string
				return fmt.Errorf("Error scanning reference contig %s: %s", name, err.Error())
			}

			refLength += len(ref)

			// calls is a chunk of calls that can be min(remaining sample contig length, defaultBufSize)
			sampleChunk, err = analyses.ReadPositions(len(ref))
			if err != nil {
				// TODO return a type that implements the error interface instead of a raw error string
				return fmt.Errorf("Error reading contig %s: %s", name, err.Error())
			}

			chunk := Chunk{
				contigName:    name,
				positionsChan: make(chan []*Position, 1),
			}

			queue <- chunk

			go analyzePositions(chunk.positionsChan, statsChan, ref, dup, sampleChunk, minCoverage, minProportion)
		}
	}

	return err
}

// TODO: rename ch to positionsChan
// Get 1 statsPool, 1 positionsPool, maxBufSize positionPool
// Put 1 callsPool
func analyzePositions(ch chan []*Position, statsChan chan SampleStats, ref, dup []byte, sampleChunks *SampleChunk, minCoverage int, minProportion float64) {
	defer close(ch)
	/*
		defer func() {
			fmt.Println("Shutdown NumGoroutine", runtime.NumGoroutine())
			close(ch)
		}()
	*/

	stats := statsPool.Get().(SampleStats)
	// Clear old values if this is a recycled object.
	stats.Reset()
	// TODO: stop reslicing?
	positions := positionsPool.Get().([]*Position)[:len(ref)]
	pattern := NewPattern()

	for i := range ref {
		refCall := ToUpper(ref[i])
		isReferenceDuplicated := len(dup) != 0 && dup[i] == '1'

		positions[i] = sampleChunks.Compare(stats, pattern, refCall, i, isReferenceDuplicated, minCoverage, minProportion)
	}
	callsPool.Put(sampleChunks.calls)
	ch <- positions
	statsChan <- stats
}

type SampleAnalysis interface {
	Name() string
	SeekContig(name string) error
	ReadPositions(n int) ([]byte, []int, []float64, error)
}

// SampleAnalyses implements sort.Interface for []SampleAnalysis based on
// the identifier field.
type SampleAnalyses []SampleAnalysis

// Get 1 callsPool
// TODO: struct { calls, cov, prop }
func (s SampleAnalyses) ReadPositions(n int) (*SampleChunk, error) {
	sampleChunk := sampleChunkPool.Get().(*SampleChunk)

	for i, analysis := range s {
		calls, coverages, proportions, err := analysis.ReadPositions(n)
		// FIXME: EOF will be common, but silently ignoring it might not be good either.
		if err != nil {
			log.Fatal(err)
		}
		//switch err {
		//default:
		//log.Fatal(err)
		//sampleChunkPool.Put(sampleChunk)
		// TODO: Panic?
		//return nil, err
		//case bufio.ErrBufferFull:
		//	log.Fatal(err)
		//	break
		//case nil, io.EOF:
		sampleChunk.calls[i] = sampleChunk.calls[i][:len(calls)]
		copy(sampleChunk.calls[i], calls)
		sampleChunk.coverages[i] = sampleChunk.coverages[i][:len(coverages)]
		copy(sampleChunk.coverages[i], coverages)
		sampleChunk.proportions[i] = sampleChunk.proportions[i][:len(proportions)]
		copy(sampleChunk.proportions[i], proportions)
		//}
	}
	return sampleChunk, nil
}

/**
 * SeekContig moves the io.Reader for each analysis onto the first position
 * of the contig with the given name, or to io.EOF if it does not exist.
 */
func (s SampleAnalyses) SeekContig(name string) (err error) {
	for _, analysis := range s {
		if err = analysis.SeekContig(name); err != nil {
			return err
		}
	}
	return nil
}

func (s SampleAnalyses) Len() int {
	return len(s)
}

func (s SampleAnalyses) Less(i, j int) bool {
	return s[i].Name() < s[j].Name()
}

func (s SampleAnalyses) Swap(i, j int) {
	s[i], s[j] = s[j], s[i]
}

func NewSampleAnalyses(files ...NaspFile) SampleAnalyses {
	analyses := make(SampleAnalyses, len(files))

	log.Printf("Indexing %d files\n", len(files))

	ch := make(chan SampleAnalysis)
	defer close(ch)

	for i := range files {
		go func(ch chan SampleAnalysis, naspFile NaspFile) {
			var analysis SampleAnalysis
			var err error

			switch {
			default:
				log.Fatalf("Unknown sample analysis: '%s' expecting *.fasta or *.vcf\n", naspFile.Filepath)
			case strings.HasSuffix(naspFile.Filepath, "fasta"):
				// NOTE: This will accept any path matching *fasta$
				// ex: .fasta, .frankenfasta, notReallyAfasta
				if analysis, err = NewFasta(naspFile, true); err != nil {
					log.Fatal(err)
				}
			case strings.HasSuffix(naspFile.Filepath, "vcf"):
				// TODO: Error handling
				if analysis, err = NewVcf(naspFile); err != nil {
					log.Fatal(err)
				}
			}

			log.Printf("\t%s %s\n", naspFile.Filepath, time.Now().Sub(t0))
			ch <- analysis
		}(ch, files[i])
	}

	for i := range analyses {
		analyses[i] = <-ch
	}

	return analyses
}

type SampleChunk struct {
	calls       [][]byte
	coverages   [][]int
	proportions [][]float64
}

func (s *SampleChunk) Compare(stats SampleStats, pattern *Pattern, refCall byte, refCallIdx int, isDup bool, minCoverage int, minProportion float64) *Position {
	var anyCalled, anyCoverage, anyProportion, anyQuality, anyRef, anySnp, anyDegen bool
	allCalled, allCoverage, allProportion, allQuality, allRef, allSnp, allDegen := true, true, true, true, true, true, true
	var call byte

	pattern.Reset()

	position := positionPool.Get().(*Position)
	position.Reset()

	switch refCall {
	default:
	case 'G', 'A', 'T', 'C':
		position.isReferenceGATC = true
	}
	position.callStr[0] = refCall
	position.maskedCallStr[0] = refCall
	position.pattern[0] = pattern.Next(refCall, true)

	position.isReferenceDuplicated = isDup

	for j := range s.calls {
		isPassCoverage, isPassProportion := false, false

		if len(s.calls[j]) < refCallIdx+1 {
			// Fasta
			call = 'X'
		} else {
			call = ToUpper(s.calls[j][refCallIdx])
		}
		position.callStr[j+1] = call

		if len(s.coverages[j]) == 0 {
			// Fasta
			isPassCoverage = true
			position.passedCoverage++
			position.passedDepthFilter[j] = '-'
		} else if call == 'X' {
			position.passedDepthFilter[j] = '?'
			position.isAllPassCoverage = false
		} else if s.coverages[j][refCallIdx] == -1 {
			isPassCoverage = true
			position.passedCoverage++
			position.passedDepthFilter[j] = '-'
		} else if s.coverages[j][refCallIdx] >= minCoverage {
			isPassCoverage = true
			position.passedCoverage++
			position.passedDepthFilter[j] = 'Y'
		} else {
			position.isAllPassCoverage = false
			position.passedDepthFilter[j] = 'N'
		}

		if len(s.proportions[j]) == 0 {
			// Fasta
			isPassProportion = true
			position.passedProportion++
			position.passedProportionFilter[j] = '-'
		} else if call == 'X' {
			// Missing VCF position.
			position.passedProportionFilter[j] = '?'
			position.isAllPassProportion = false
		} else if s.proportions[j][refCallIdx] == -1 {
			// Cannot determine due to no value specified.
			isPassProportion = true
			position.passedProportion++
			position.passedProportionFilter[j] = '-'
		} else if s.proportions[j][refCallIdx] >= minProportion {
			// Pass
			isPassProportion = true
			position.passedProportion++
			position.passedProportionFilter[j] = 'Y'
		} else {
			// Fail
			position.passedProportionFilter[j] = 'N'
			position.isAllPassProportion = false
		}

		if !(isPassCoverage && isPassProportion) {
			// TODO: Sample consensus
			// We can short-circuit the consensus check here because it has already failed
			position.isAllPassConsensus = false
			//position.isAllQualityBreadth = false
		}

		isDegen, wasCalled := position.CheckCall(j, call, isPassProportion && isPassCoverage)

		// TODO: extract if/else and isQualty, isRef, isSnp bools as function
		isQuality := true // assume true until proven otherwise
		isRef, isSnp := false, false
		if wasCalled && position.isReferenceGATC && isPassProportion && isPassCoverage {
			if position.isReferenceDuplicated {
				isQuality = false
				position.isAllQualityBreadth = false
			}

			if isDegen {
				position.calledDegen++
			} else if refCall == call {
				isRef = true
				position.calledReference++
			} else {
				isSnp = true
				position.calledSnp++
			}

		} else {
			isQuality = false
			position.isAllQualityBreadth = false
		}

		if !position.isReferenceDuplicated {
			// 2015-08-18 coverage and proportion sample stats do not include
			// duplicate regions to be consistent with the other stats.
			if isPassCoverage {
				stats[j].passedCoverageFilter++
				anyCoverage = true
			} else {
				allCoverage = false
			}

			if isPassProportion {
				stats[j].passedProportionFilter++
				anyProportion = true
			} else {
				allProportion = false
			}

			if wasCalled {
				stats[j].wasCalled++
				anyCalled = true
			} else {
				allCalled = false
			}

			if isQuality {
				stats[j].qualityBreadth++
				anyQuality = true
			} else {
				allQuality = false
			}

			if isRef {
				stats[j].calledReference++
				anyRef = true
				allSnp, allDegen = false, false
			} else if isSnp {
				// Called A/C/G/T and doesn't match the reference
				stats[j].calledSnp++
				position.isAnySnp = true
				position.isMissingMatrix = true
				anySnp = true
				allRef, allDegen = false, false
			} else if isDegen {
				stats[j].calledDegen++
				anyDegen = true
				allRef, allSnp = false, false
			}
		} else {
			allCoverage, allProportion, allCalled, allQuality, allRef, allSnp, allDegen = false, false, false, false, false, false, false
		}

		// *.tsv Pattern Column
		position.pattern[j+1] = pattern.Next(call, isPassCoverage && isPassProportion)

		position.isBestSnp = position.isAllQualityBreadth && position.calledSnp > 0
	}

	// Whole Genome Any/All Sample Stats summary are the last two elements of
	// of the stats array.
	stats[len(stats)-2].Increment(anyCalled, anyCoverage, anyProportion, anyQuality, anyRef, anySnp, anyDegen)
	stats[len(stats)-1].Increment(allCalled, allCoverage, allProportion, allQuality, allRef, allSnp, allDegen)

	return position
}

/**
 * ToUpper is reduced from the standard library version which includes unicode
 * support and consequently additional memory allocation.
 * It promotes an ASCII byte to uppercase.
 */
func ToUpper(b byte) byte {
	if 'a' <= b && b <= 'z' {
		b -= 'a' - 'A'
	}
	return b
}
