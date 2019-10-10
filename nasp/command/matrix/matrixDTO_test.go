package matrix

import (
	"encoding/xml"
	"reflect"
	"testing"
)

const data = `
<?xml version="1.0" ?>
<matrix_data>
    <parameters>
        <minimum-coverage>10</minimum-coverage>
        <reference-fasta>/home/jtravis/NASP/examples/example_1/nasp_results/reference/reference.fasta</reference-fasta>
        <matrix-folder>/home/jtravis/NASP/examples/example_1/nasp_results/matrices</matrix-folder>
        <minimum-proportion>0.9</minimum-proportion>
        <stats-folder>/home/jtravis/NASP/examples/example_1/nasp_results/statistics</stats-folder>
        <reference-dups>/home/jtravis/NASP/examples/example_1/nasp_results/reference/duplicates.txt</reference-dups>
    </parameters>
    <files>
        <frankenfasta aligner="nucmer" name="example_1">/home/jtravis/NASP/examples/example_1/nasp_results/external/example_1.frankenfasta</frankenfasta>
        <vcf aligner="Bowtie2" name="example_1_L001-bowtie2-gatk" snpcaller="GATK">/home/jtravis/NASP/examples/example_1/nasp_results/gatk/example_1_L001-bowtie2-gatk.vcf</vcf>
    </files>
</matrix_data>
`

func TestParseExampleMatrixDto(t *testing.T) {
	var actual Dto
	if err := xml.Unmarshal([]byte(data), &actual); err != nil {
		t.Fatal(err)
	}

	expect := Dto{
		Parameters: Parameters{
			MinCoverage:   10,
			MinProportion: 0.9,
			MatrixFolder:  "/home/jtravis/NASP/examples/example_1/nasp_results/matrices",
			StatsFolder:   "/home/jtravis/NASP/examples/example_1/nasp_results/statistics",
			Reference:     "/home/jtravis/NASP/examples/example_1/nasp_results/reference/reference.fasta",
			Duplicates:    "/home/jtravis/NASP/examples/example_1/nasp_results/reference/duplicates.txt",
		},
		Fastas: []*NaspFile{
			{Name: "example_1", Aligner: "nucmer", Snpcaller: "", Filepath: "/home/jtravis/NASP/examples/example_1/nasp_results/external/example_1.frankenfasta"},
		},
		Vcfs: []*NaspFile{
			{Name: "example_1_L001-bowtie2-gatk", Aligner: "Bowtie2", Snpcaller: "GATK", Filepath: "/home/jtravis/NASP/examples/example_1/nasp_results/gatk/example_1_L001-bowtie2-gatk.vcf"},
		},
	}

	if !reflect.DeepEqual(expect, actual) {
		t.Fatalf("expect:\n%#v\nactual:\n%#v\n", expect, actual)
	}
}

func TestNewDto(t *testing.T) {
	t.Log("TODO: incomplete")
	var tests = []struct {
		dtoPath       string
		refPath       string
		dupPath       string
		matrixPath    string
		statsPath     string
		minCoverage   int
		minProportion float64
		files         []string
	}{
		{},
	}

	for _, tt := range tests {
		dto, err := NewDto(tt.dtoPath, tt.refPath, tt.dupPath, tt.matrixPath, tt.statsPath, tt.minCoverage, tt.minProportion, tt.files)
		if err != nil {
			t.Fatal(err)
		}
		t.Logf("%#v\n", dto)
	}
}

func TestNewNaspFileFromFilepath(t *testing.T) {
	var tests = []struct {
		input  string
		expect *NaspFile
	}{
		{
			input:  "",
			expect: &NaspFile{},
		}, {
			input: "example_1_L001.frankenfasta",
			expect: &NaspFile{
				Name:      "example_1_L001",
				Aligner:   "nucmer",
				Snpcaller: "",
				Filepath:  "example_1_L001.frankenfasta",
			},
		}, {
			input: "example_1_L001-bowtie2-gatk.vcf",
			expect: &NaspFile{
				Name:      "example_1_L001",
				Aligner:   "Bowtie2",
				Snpcaller: "GATK",
				Filepath:  "example_1_L001-bowtie2-gatk.vcf",
			},
		}, {
			input: "example_1_L001.vcf",
			expect: &NaspFile{
				Name:      "example_1_L001",
				Aligner:   "pre-aligned",
				Snpcaller: "pre-called",
				Filepath:  "example_1_L001.vcf",
			},
		},
	}

	for _, tt := range tests {
		actual := NewNaspFileFromFilepath(tt.input)
		if !reflect.DeepEqual(tt.expect, actual) {
			t.Errorf("expect: %#v actual: %#v", tt.expect, actual)
		}
	}
}
