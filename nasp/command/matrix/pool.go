package matrix

import "sync"

var callsPool = sync.Pool{
	New: func() interface{} {
		// Allocate 1 array instead of N. Unproven assumption that a single
		// contiguous chunk of memory is more efficient.
		b := make([]byte, NUM_SAMPLES*defaultBufSize)
		calls := make([][]byte, NUM_SAMPLES)
		for i := range calls {
			calls[i] = b[i*defaultBufSize : (i+1)*defaultBufSize]
		}
		return calls
	},
}

var positionsPool = sync.Pool{
	New: func() interface{} {
		return make([]*Position, defaultBufSize)
	},
}

type Position struct {
	// General Stats
	isAllCalled           bool
	isReferenceGATC       bool
	isReferenceDuplicated bool
	isAllPassCoverage     bool
	isAllPassProportion   bool
	isAllPassConsensus    bool
	isAllQualityBreadth   bool
	isAnySnp              bool
	isBestSnp             bool
	isMissingMatrix       bool

	//all_sample_stats=all_sample_stats,

	// Missing Data Matrix condition - at least one SampleAnalysis passes quality_breadth and is a SNP.
	//is_missing_matrix=is_missing_matrix,

	// NASP Master Matrix
	// Counters
	wasCalled        int64
	calledReference  int64
	calledSnp        int64
	calledDegen      int64
	passedCoverage   int64
	passedProportion int64
	a                int64
	c                int64
	g                int64
	t                int64
	n                int64

	// Strings
	callStr                []byte
	maskedCallStr          []byte
	callWasMade            []byte
	passedDepthFilter      []byte
	passedProportionFilter []byte
	pattern                []byte
}

func (p *Position) Reset() {
	p.isAllCalled = true
	p.isReferenceGATC = false
	p.isReferenceDuplicated = false
	p.isAllPassCoverage = true
	p.isAllPassProportion = true
	p.isAllPassConsensus = true
	p.isAllQualityBreadth = true
	p.isAnySnp = false
	p.isBestSnp = true
	p.isMissingMatrix = false

	p.wasCalled = 0
	p.calledReference = 0
	p.calledSnp = 0
	p.calledDegen = 0
	p.passedCoverage = 0
	p.passedProportion = 0
	p.a = 0
	p.c = 0
	p.g = 0
	p.t = 0
	p.n = 0
}

func (position *Position) CheckCall(sampleIdx int, call byte, isPassFilters bool) (bool, bool) {
	// It can be called anything so long as it was called something
	// X and N
	isDegen := false
	wasCalled := true

	switch call {
	case 'G':
		position.g++
	case 'A':
		position.a++
	case 'T':
		position.t++
	case 'C':
		position.c++
	case 'X', 'N':
		wasCalled = false
		fallthrough
	default:
		position.n++
		isDegen = true
		position.isAllPassConsensus = false
		position.isAllQualityBreadth = false
	}

	if !wasCalled || isPassFilters && !isDegen {
		position.maskedCallStr[sampleIdx+1] = call
	} else {
		position.maskedCallStr[sampleIdx+1] = 'N'
		position.isAllPassConsensus = false
		position.isAllQualityBreadth = false
	}

	if wasCalled {
		position.wasCalled++
		position.callWasMade[sampleIdx] = 'Y'
	} else {
		position.callWasMade[sampleIdx] = 'N'
		position.isAllCalled = false
		position.isAllQualityBreadth = false
	}

	return isDegen, wasCalled
}

var positionPool = sync.Pool{
	New: func() interface{} {
		// Allocate 1 array instead of N. Unproven assumption that a single
		// contiguous chunk of memory is more efficient.
		b := make([]byte, 6*NUM_SAMPLES+3)
		return &Position{
			isAllCalled: true,
			//isReferenceGATC:      true,
			//isReferenceDuplicated: true,
			isAllPassCoverage:   true,
			isAllPassProportion: true,
			isAllPassConsensus:  true,
			isAllQualityBreadth: true,
			isBestSnp:           true,
			//isMissingMatrix: true,

			callWasMade:            b[:NUM_SAMPLES],
			passedDepthFilter:      b[NUM_SAMPLES : 2*NUM_SAMPLES],
			passedProportionFilter: b[2*NUM_SAMPLES : 3*NUM_SAMPLES],
			callStr:                b[3*NUM_SAMPLES : 4*NUM_SAMPLES+1],
			maskedCallStr:          b[4*NUM_SAMPLES+1 : 5*NUM_SAMPLES+2],
			pattern:                b[5*NUM_SAMPLES+2:],
		}
	},
}

var statsPool = sync.Pool{
	New: func() interface{} {
		// stats has an element for each sample +2 for the Whole Genome any/all
		// summary as the last two elements.
		return make(SampleStats, NUM_SAMPLES+2)
	},
}

var sampleChunkPool = sync.Pool{
	New: func() interface{} {
		// Allocate 1 array instead of N. Unproven assumption that a single
		// contiguous chunk of memory is more efficient.
		calls := make([]byte, NUM_SAMPLES*defaultBufSize)
		coverages := make([]int, NUM_SAMPLES*defaultBufSize)
		proportions := make([]float64, NUM_SAMPLES*defaultBufSize)

		sampleChunk := &SampleChunk{
			calls:       make([][]byte, NUM_SAMPLES),
			coverages:   make([][]int, NUM_SAMPLES),
			proportions: make([][]float64, NUM_SAMPLES),
		}

		for i := 0; i < NUM_SAMPLES; i++ {
			sampleChunk.calls[i] = calls[i*defaultBufSize : (i+1)*defaultBufSize]
			sampleChunk.coverages[i] = coverages[i*defaultBufSize : (i+1)*defaultBufSize]
			sampleChunk.proportions[i] = proportions[i*defaultBufSize : (i+1)*defaultBufSize]
		}

		return sampleChunk
	},
}
