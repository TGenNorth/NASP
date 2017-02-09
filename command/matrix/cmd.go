// TODO: LICENSE
// Jason Travis
package matrix

import (
	"runtime"

	"github.com/TGenNorth/NASP/command"
)

var cmd = &command.Command{
	UsageLine: "matrix --dto-file matrix_dto.xml",
	Short:     "aggregate sample calls into matrices and stats",
	Long: `
Matrix aggregates the VCF and Frankenfasta files from the NASP pipeline into .tsv matrix and stats file summaries
which may be parsed for further analysis.

Given the --dto-file flag, all other flags are optional overrides.
A reference fasta and vcf/frankenfasta files are always required, either from
the dto or listed on the command line.

--dto-file           Path to the matrix_dto.xml file
--num-threads        Max number of CPUs that can be executing simultaneously (default: all)
--reference-fasta    Path to the reference.fasta against which samples are compared
--reference-dups     Path to the duplicates.txt file marking duplicated positions
--stats-folder       Path to the output folder for statistics (default: ./)
--matrix-folder      Path to the output folder for matrices (default: ./)
--minimum-coverage   Filter positions below this coverage/depth threshold (default: 0)
--minimum-proportion Filter positions below this proportion threshold (default: 0.0)
--withallrefpos      Include the withallrefpos.tsv matrix
	`,
}

var (
	numThreads    int
	dtoFile       string
	refPath       string
	dupPath       string
	statsPath     string
	matrixPath    string
	minCoverage   int
	minProportion float64
	withAllRefPos bool
)

func init() {
	// break init cycle
	cmd.Run = func(cmd *command.Command, args []string) error {
		return Run(numThreads, dtoFile, refPath, dupPath, statsPath, matrixPath, minCoverage, minProportion, args...)
	}
	cmd.Flag.IntVar(&numThreads, "num-threads", runtime.NumCPU(), "Max number of CPUs that can be executing simultaneously")
	cmd.Flag.StringVar(&dtoFile, "dto-file", "", "Path to the matrix_dto.xml file")
	cmd.Flag.StringVar(&refPath, "reference-fasta", "", "Path to the reference.fasta against which samples are compared")
	cmd.Flag.StringVar(&dupPath, "reference-dups", "", "Path to the duplicates.txt file marking duplicated positions")
	cmd.Flag.StringVar(&statsPath, "stats-folder", "", "Path to the output folder for statistics")
	cmd.Flag.StringVar(&matrixPath, "matrix-folder", "", "Path to the output folder for matrices")
	cmd.Flag.IntVar(&minCoverage, "minimum-coverage", -1, "Filter positions below this coverage/depth threshold")
	cmd.Flag.Float64Var(&minProportion, "minimum-proportion", -1.0, "Filter positions below this proportion threshold")
	cmd.Flag.BoolVar(&withAllRefPos, "withallrefpos", false, "Include the withallrefpos matrix")

	command.Register(cmd)
}
