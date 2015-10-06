package matrix

import (
	"encoding/xml"
	"io/ioutil"
	"log"
	"os"
	"path/filepath"
	"strings"
)

type Parameters struct {
	MinCoverage        int     `xml:"parameters>minimum-coverage"`
	MinProportion      float64 `xml:"parameters>minimum-proportion"`
	Reference          string  `xml:"parameters>reference-fasta"`
	Duplicates         string  `xml:"parameters>reference-dups"`
	MatrixFolder       string  `xml:"parameters>matrix-folder"`
	StatsFolder        string  `xml:"parameters>stats-folder"`
	FilterMatrixFormat string  `xml:"parameters>filter-matrix-format"`
}

type NaspFile struct {
	Name      string `xml:"name,attr"`
	Aligner   string `xml:"aligner,attr"`
	Snpcaller string `xml:"snpcaller,attr"`
	Filepath  string `xml:",chardata"`
}

func NewNaspFileFromFilepath(path string) NaspFile {
	var aligner, snpcaller string

	name := filepath.Base(path)
	extension := filepath.Ext(name)

	// Trim extension
	name = name[0 : len(name)-len(extension)]

	switch extension {
	default:
	case ".frankenfasta":
		aligner = "nucmer"
	case ".vcf":
		// Trim SnpCaller
		switch {
		default:
			snpcaller = "pre-called"
		case strings.HasSuffix(name, "-gatk"):
			snpcaller = "GATK"
			name = name[:len(name)-len("-gatk")]
		case strings.HasSuffix(name, "-varscan"):
			snpcaller = "VarScan"
			name = name[:len(name)-len("-varscan")]
		case strings.HasSuffix(name, "-samtools"):
			snpcaller = "SAMtools"
			name = name[:len(name)-len("-samtools")]
		case strings.HasSuffix(name, "-solsnp"):
			snpcaller = "SolSNP"
			name = name[:len(name)-len("-solsnp")]
		}

		// Trim Aligner
		switch {
		default:
			aligner = "pre-aligned"
		case strings.HasSuffix(name, "-novo"):
			aligner = "Novoalign"
			name = name[:len(name)-len("-novo")]
		case strings.HasSuffix(name, "-bowtie2"):
			aligner = "Bowtie2"
			name = name[:len(name)-len("-bowtie2")]
		case strings.HasSuffix(name, "-bwamem"):
			aligner = "BWA-mem"
			name = name[:len(name)-len("-bwamem")]
		//case strings.HasSuffix(name, "-bwamem"):
		//	// TODO: get sampe suffix
		//	aligner = "BWA-sampe"
		//	name = name[:len(name)-len("-bwasampe")]
		case strings.HasSuffix(name, "-snap"):
			aligner = "SNAP"
			name = name[:len(name)-len("-snap")]
		}
	}

	return NaspFile{
		Aligner:   aligner,
		Filepath:  path,
		Name:      name,
		Snpcaller: snpcaller,
	}
}

type Dto struct {
	Parameters `xml:"parameters"`
	Fastas     []NaspFile `xml:"files>frankenfasta"`
	Vcfs       []NaspFile `xml:"files>vcf"`
	// All files, including files listed on the commandline
	AllFiles []NaspFile
}

// - Commandline arguments supersede dto file values.
// - minCoverage and minProportion must be -1 to indicate
func NewDto(dtoPath, refPath, dupPath, matrixPath, statsPath string, minCoverage int, minProportion float64, files []string) (*Dto, error) {
	dto := &Dto{}

	// Parse the matrix_dto.xml
	if dtoPath != "" {
		file, err := os.Open(dtoPath)
		if err != nil {
			return nil, err
		}
		defer file.Close()

		data, err := ioutil.ReadAll(file)
		if err != nil {
			return nil, err
		}

		if err := xml.Unmarshal(data, dto); err != nil {
			return nil, err
		}
	}

	// Commandline arguments supersede matrix_dto.xml values.
	if refPath != "" {
		dto.Reference = refPath
	}

	if dupPath != "" {
		dto.Duplicates = dupPath
	}

	if statsPath != "" {
		dto.StatsFolder = statsPath
	}

	if matrixPath != "" {
		dto.MatrixFolder = matrixPath
	}

	if minCoverage != -1 {
		dto.MinCoverage = minCoverage
	}

	if minProportion != -1 {
		dto.MinProportion = minProportion
	}

	// TODO: what if the stats/matrix folder is undefined?

	if dto.StatsFolder != "" {
		if err := os.MkdirAll(dto.StatsFolder, 0777); err != nil {
			return nil, err
		}
	}

	if dto.MatrixFolder != "" {
		if err := os.MkdirAll(dto.MatrixFolder, 0777); err != nil {
			return nil, err
		}
	}

	// The reference file is required and must exist.
	if _, err := os.Stat(dto.Reference); err != nil {
		if os.IsNotExist(err) {
			return nil, errRequired("a reference fasta")
		}
		return nil, err
	}

	// The duplicates file is optional. If specified, make sure it exists.
	if _, err := os.Stat(dto.Duplicates); dto.Duplicates != "" && err != nil {
		return nil, errNotExist(dto.Duplicates)
	}

	if err := dto.mergeAllFiles(files); err != nil {
		return nil, err
	}

	return dto, nil
}

func (d *Dto) mergeAllFiles(filepaths []string) error {
	// Ignore duplicate filepaths
	// While this is not the same as ensuring no two files are the same, it is
	// a low cost check. Examples that would fool this check:
	// Two symbolic links to the same file or two copies of the same file with
	// different filenames or different locations.
	uniqueFilepaths := make(map[string]NaspFile)

	// TODO: avoid allocation if all files are VCFs or Fastas
	for i := range filepaths {
		if _, err := os.Stat(filepaths[i]); os.IsNotExist(err) {
			return errNotExist(filepaths[i])
		}
		if _, ok := uniqueFilepaths[filepaths[i]]; ok {
			log.Printf("Duplicate filepath: '%s'", filepaths[i])
		} else {
			uniqueFilepaths[filepaths[i]] = NewNaspFileFromFilepath(filepaths[i])
		}
	}
	for i := range d.Vcfs {
		if _, err := os.Stat(d.Vcfs[i].Filepath); os.IsNotExist(err) {
			return errNotExist(d.Vcfs[i].Filepath)
		}
		if _, ok := uniqueFilepaths[d.Vcfs[i].Filepath]; ok {
			log.Printf("Duplicate filepath: '%s'", d.Vcfs[i].Filepath)
		} else {
			uniqueFilepaths[d.Vcfs[i].Filepath] = d.Vcfs[i]
		}
	}
	for i := range d.Fastas {
		if _, err := os.Stat(d.Fastas[i].Filepath); os.IsNotExist(err) {
			return errNotExist(d.Fastas[i].Filepath)
		}
		if _, ok := uniqueFilepaths[d.Fastas[i].Filepath]; ok {
			log.Printf("Duplicate filepath: '%s'", d.Fastas[i].Filepath)
		} else {
			uniqueFilepaths[d.Fastas[i].Filepath] = d.Fastas[i]
		}
	}

	d.AllFiles = make([]NaspFile, len(uniqueFilepaths))
	i := 0
	for _, v := range uniqueFilepaths {
		d.AllFiles[i] = v
		i++
	}

	// The program will panic with a negative workgroup error if there is not
	// at least one sample analysis
	if len(uniqueFilepaths) == 0 {
		return errRequired("at least one sample analysis (vcf/frankenfasta)")
	}

	return nil
}

type errNotExist string

func (e errNotExist) Error() string {
	return "nasp: file does not exist '" + string(e) + "'"
}

type errRequired string

func (e errRequired) Error() string {
	return "nasp: " + string(e) + " must be specified either in the matrix_dto.xml or on the command line"
}
