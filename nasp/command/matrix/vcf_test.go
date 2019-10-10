package matrix

import (
	"bytes"
	"fmt"
	"log"
	"path/filepath"
	"runtime"
	"testing"
)

func init() {
	var err error
	paths, err = filepath.Glob("./klebs/gatk/*.vcf")
	if err != nil {
		log.Fatal(err)
	}
	runtime.GOMAXPROCS(runtime.NumCPU())
}

var paths []string

/*
func BenchmarkAllVcfs(b *testing.B) {
	analyses := make([]*Vcf, len(paths))

	t0 := time.Now()

	for i := range paths {
		analyses[i] = NewVcf("", "", paths[i])
		fmt.Println(time.Now().Sub(t0))
	}
}

func BenchmarkAllVcfsGoroutineUnBuffered(b *testing.B) {
	analyses := make([]*Vcf, len(paths))

	t0 := time.Now()

	ch := make(chan *Vcf)
	for i := range paths {
		go func(ch chan *Vcf, path string) {
			ch <- NewVcf("", "", path)
			fmt.Println(time.Now().Sub(t0))
		}(ch, paths[i])
	}

	for i := range analyses {
		analyses[i] = <-ch
	}
}

func BenchmarkAllVcfsGoroutineBuffered(b *testing.B) {
	analyses := make([]*Vcf, len(paths))

	t0 := time.Now()

	ch := make(chan *Vcf, 1)
	for i := range paths {
		go func(ch chan *Vcf, path string) {
			ch <- NewVcf("", "", path)
			fmt.Println(time.Now().Sub(t0))
		}(ch, paths[i])
	}

	for i := range analyses {
		analyses[i] = <-ch
	}
}

func BenchmarkAllVcfsGoroutineChunks(b *testing.B) {
	fmt.Println(runtime.GOMAXPROCS(0))
	t0 := time.Now()

	analyses := make([]*Vcf, len(paths))

	fmt.Printf("Indexing %d files\n", len(paths))

	chunkSize := len(paths) / 4

	fmt.Printf("Breaking into %d size chunks\n", chunkSize)

	ch := make(chan []*Vcf, 4)
	for i := 0; i+chunkSize < len(paths); i += chunkSize {
		go func(ch chan []*Vcf, p []string) {
			chunk := make([]*Vcf, len(p))
			for j := range p {
				chunk[j] = NewVcf("", "", p[j])
			}

			ch <- chunk
			fmt.Println(len(p), time.Now().Sub(t0))
		}(ch, paths[i:i+chunkSize])
	}

	for i := 0; i < 4; i++ {
		c := <-ch
		for j, v := range c {
			analyses[i*chunkSize+j] = v
		}
	}
}

func BenchmarkNewVcfSingle(b *testing.B) {
	NewVcf("", "", paths[0])
}

func BenchmarkVcfNewSampleAnalyses(b *testing.B) {
	NewSampleAnalyses(paths...)
}

func BenchmarkVcfReadPositions(b *testing.B) {
	v, err := NewVcf("", "", "/Users/jtravis/workspace/src/github.com/corburn/gonasp/klebs/gatk/Kpneumo-Sp29_CTTGTA_L007-novo-gatk.vcf")
	if err != nil {
		b.Fatal(err)
	}
	if err := v.SeekContig("KP1_0026"); err != nil {
		b.Fatal(err)
	}

	var calls []byte
	var coverages []int
	var proportions []float64
	var err error
	for ; err == nil; calls, coverages, proportions, err = v.ReadPositions(defaultBufSize) {
		fmt.Printf("%s, %v, %v", calls, coverages, proportions)
	}
	fmt.Println(err)
	//fmt.Println(v)
}
*/

const (
	vcfData = `metadata
	ContigA contains a gap before the first record
	ContigB contains a gap in the middle
	ContigC contains a implicit gap at the end
	ContigD contains a mix of gaps before/middle/after
	metadata
	#CHROM  POS ID  REF ALT QUAL    FILTER  INFO    FORMAT  sample
	ContigA    3   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	ContigA    4   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	ContigA    5   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	ContigB    1   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	pontigB    2   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	ContigB    5   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	ContigB    1   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	ContigB    2   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	ContigC    1   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	ContigC    2   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	ContigD    2   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	ContigD    4   .   A   .   139 .   AN=1;DP=25;MQ=70.00;MQ0=0   GT:DP:MLPSAC:MLPSAF 0:25
	`
)

func TestVcfReadPositionsBeforeSeek(t *testing.T) {
	t.Skip("TODO: Not implemented")
	v, err := NewVcf(&NaspFile{Filepath: "/Users/jtravis/workspace/src/github.com/corburn/gonasp/klebs/gatk/Kpneumo-Sp29_CTTGTA_L007-novo-gatk.vcf"})
	if err != nil {
		t.Fatal(err)
	}
	fmt.Println(v.ReadPositions(42))
}

/*
func TestVcfGapBeforeRecord(t *testing.T) {
	file, err := ioutil.TempFile("", "")
	if err != nil {
		t.Fatal(err)
	}
	defer func() {
		file.Close()
		os.Remove(file.Name())
	}()

}
*/

func TestVcfRecordSplit(t *testing.T) {
	line := []byte("KP1_0026\t12\t.\tG\t.\t117\t.\tAN=1;DP=20;MQ=70.00;MQ0=0\tGT:DP:MLPSAC:MLPSAF\t0:20")

	expect := [][]byte{
		[]byte("KP1_0026"),
		[]byte("12"),
		[]byte("."),
		[]byte("G"),
		[]byte("."),
		[]byte("117"),
		[]byte("."),
		[]byte("AN=1;DP=20;MQ=70.00;MQ0=0"),
		[]byte("GT:DP:MLPSAC:MLPSAF"),
		[]byte("0:20"),
	}

	record := make(vcfRecord, 10)
	record.Split(line)

	for i := range record {
		if !bytes.Equal(expect[i], record[i]) {
			t.Errorf("VcfRecordSplit(%s): expect %s, actual %s", line, expect[i], record[i])
		}
	}
}

func TestVcfRecordPosition(t *testing.T) {
	var tests = []struct {
		input    []byte // input
		expected int    // expected result
	}{
		// Base case.
		{[]byte("42"), 42},
	}

	record := make(vcfRecord, 10)
	for _, tt := range tests {
		record[pos] = tt.input
		actual := record.Position()
		if actual != tt.expected {
			t.Errorf("VcfRecordPosition(%s): expected %d, actual %d", tt.input, tt.expected, actual)
		}
	}
}

func TestVcfRecordCoverageInSampleColumn(t *testing.T) {
	var tests = []struct {
		format   []byte
		sample   []byte
		expected int
	}{
		// Base case.
		{[]byte(""), []byte(""), -1},
		// Ideal scenario.
		{[]byte("GT:DP"), []byte("0:42"), 42},
		// GT not present
		{[]byte("DP"), []byte("42"), 42},
		// Preceeded by at least one sub-field.
		{[]byte("FOO:DP"), []byte("24:42"), 42},
		{[]byte("FOO"), []byte("24:"), -1},
		// Another field happens to end with DP.
		{[]byte("XDP:DP"), []byte("24:42"), 42},
		// The field does not exist.
		{[]byte("FOO:XDP"), []byte("24:42"), -1},
	}

	record := make(vcfRecord, 10)
	for _, tt := range tests {
		record[format] = tt.format
		record[sample] = tt.sample
		actual := record.Coverage()
		if actual != tt.expected {
			t.Errorf("VcfRecordCoverage(%s): expected %d, actual %d\n", tt.format, tt.expected, actual)
		}
	}
}

func TestVcfRecordCoverageInInfoColumn(t *testing.T) {
	var tests = []struct {
		input    []byte // input
		expected int    // expected result
	}{
		// Base case.
		{[]byte(""), -1},
		// Ideal scenario.
		{[]byte("DP=42"), 42},
		{[]byte("ADP=42"), 42},
		// Preceeded by at least one sub-field.
		{[]byte("FOO;DP=42"), 42},
		{[]byte("FOO;ADP=42"), 42},
		{[]byte("FOO"), -1},
		// Another field happens to end with DP.
		{[]byte("XDP=24;DP=42"), 42},
		{[]byte("XADP=24;ADP=42"), 42},
		// The field does not exist.
		{[]byte("FOO=24;XDP=42"), -1},
		{[]byte("FOO=24;XADP=42"), -1},
		// It stops searching after the field is found
		{[]byte("AC=1;AF=1.00;AN=1;DP=15;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=64.30;MQ0=0;QD=11.04"), 15},
	}

	record := make(vcfRecord, 10)
	for _, tt := range tests {
		record[info] = tt.input
		actual := record.Coverage()
		if actual != tt.expected {
			t.Errorf("VcfRecordCoverage(%s): expected %d, actual %d\n", tt.input, tt.expected, actual)
		}
	}
}

func TestVcfRecordSamtoolsProportion(t *testing.T) {
	var tests = []struct {
		input    []byte
		expected float64
	}{
		// bowtie2-samtools
		/*{
			input:    []byte("500WT1_test\t1\t.\tA\tG\t5.46\t.\tDP=961;VDB=0.000000e+00;RPB=2.642707e+01;AF1=0.4999;AC1=1;DP4=80,4,421,0;MQ=27;FQ=7.8;PV4=1,1,7.4e-35,1\tGT:PL:DP:GQ\t0/1:34,0,153:2:34"),
			expected: 42.0,
		},*/
		{
			input:    []byte("500WT1_test\t1373\t.\tC\tG\t222\t.\tDP=31982;VDB=1.142963e-01;RPB=5.297841e-01;AF1=1;AC1=2;DP4=2,2,12589,13869;MQ=42;FQ=-282;PV4=1,1,1,0.11\tGT:PL:DP:GQ\t1/1:255,255,0:26462:99"),
			expected: 42.0,
		},
	}

	record := make(vcfRecord, 10)
	for _, tt := range tests {
		record.Split(tt.input)
		call := record.Call()
		coverage := record.Coverage()
		proportion := record.Proportion(call, coverage)
		if proportion != tt.expected {
			t.Errorf("test %s call: %c coverage: %d proportion: %f expected: %f\n", tt.input, call, coverage, proportion, tt.expected)
		}
	}
	//500WT1_test     1       .       A       G       5.46    .       DP=961;VDB=0.000000e+00;RPB=2.642707e+01;AF1=0.4999;AC1=1;DP4=522,0,421,0;MQ=27;FQ=7.8;PV4=1,1,7.4e-35,1        GT:PL:DP:GQ     0/1:34,0,153:943:34
}

func TestVcfRecordGatkVarscanProportion(t *testing.T) {
	var tests = []struct {
		call     byte
		ref      []byte
		format   []byte
		sample   []byte
		expected float64
	}{
		// Base case.
		//{'.', []byte("."), []byte(""), []byte(""), -1},
		// GATK
		{'.', []byte("."), []byte("GT:AD"), []byte("0:84"), 42.0},
		//{'.', []byte("."), []byte("GT:AD"), []byte("0:42,84"), 21.0},
		//{'.', []byte("."), []byte("GT:AD"), []byte("1:0,84"), 42.0},
		// It should stop searching after finding the field.
		// This position fails a .9 proportion filter if the sample column is
		// parsed first, but passes if the info column is parsed first.
		// Why is the DP value different?
		// Saureus-FMC1_L005-novo-gatk.vcf
		// NC_007790.1     671     .       CGAATA  C       827.97  .       AC=1;AF=1.00;AN=1;DP=15;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=64.30;MQ0=0;QD=11.04     GT:AD:DP:GQ:MLPSAC:MLPSAF:PL    1:0,13:14:99:1:1.00:867,0
		//{'C', []byte("T"), []byte("GT:AD:DP:GQ:MLPSAC:MLPSAF:PL"), []byte("1:44,120:164:99:1:1.00:417,0"), 60.0},
		// Saureus-FMC8_L005-novo-gatk.vcf
		// NC_007790.1	793	.	G	.	8266	.	AN=1;DP=359;MQ=70.00;MQ0=0	GT:DP:MLPSAC:MLPSAF	0:359
		//{'G', []byte("G"), []byte("GT:DP:MLPSAC:MLPSAF"), []byte("0:359"), -1.0},
	}

	record := make(vcfRecord, 10)
	for _, tt := range tests {
		record[ref] = tt.ref
		record[format] = tt.format
		record[sample] = tt.sample
		actual := record.Proportion(tt.call, 2)
		if actual != tt.expected {
			t.Errorf("test %s %c %s %s result %f\n", tt.ref, tt.call, tt.format, tt.sample, actual)
		}
	}
}

func TestVcfRecordCall(t *testing.T) {
	var tests = []struct {
		input  []byte
		expect byte
	}{
		{input: []byte("srmB\t1212\t.\tGA\tG\t257.97\t.\tAC=1;AF=1.00;AN=1;DP=7;FS=0.000;MLEAC=1;MLEAF=1.00;MQ=70.00;MQ0=0;QD=34.80;RPA=2,1;RU=A;STR\tGT:AD:DP:GQ:MLPSAC:MLPSAF:PL\t1:0,7:7:99:1:1.00:297,0"), expect: 'G'},
		{input: []byte("500WT1_test\t320\t.\tGTT\tTT\t.\tPASS\tADP=16608;WT=0;HET=0;HOM=1;NC=0\tGT:GQ:SDP:DP:RD:AD:FREQ:PVAL:RBQ:ABQ:RDF:RDR:ADF:ADR\t1/1:255:16617:16608:650:15645:94.75%:0E0:35:37:494:156:7782:7863"), expect: 'G'},
	}
	record := make(vcfRecord, 10)

	for _, tt := range tests {
		record.Split(tt.input)
		if actual := record.Call(); actual != tt.expect {
			t.Errorf("VcfRecordCall(%s): expected %c, actual %c", tt.input, tt.expect, actual)
		}
	}
}
