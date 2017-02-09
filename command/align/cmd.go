package align

import (
	"errors"
	"fmt"
	"os"
	"strings"
	"text/template"

	"github.com/TGenNorth/NASP/command"
)

var cmd = &command.Command{
	UsageLine: "align --aligner [aligner_name[,aligner_name,...]] --reference [reference.fasta] [*.fastq]",
	Short:     "aggregate sample calls into matrices and stats",
	Long: `
Align TODO

For more about <ACTION>, see 'nasp help <ACTION>'.
	`,
}

var (
	alignersFlag   aligners
	dryRunFlag     bool
	jobManagerFlag string
	referenceFlag  string
)

func init() {
	cmd.Run = runAlign
	cmd.Flag.Var(&alignersFlag, "aligner", "comma-separated list of aligners")
	cmd.Flag.BoolVar(&dryRunFlag, "dry-run", true, "write the job script, but do not run it")
	cmd.Flag.StringVar(&referenceFlag, "reference", "", "path to the reference fasta")
	cmd.Flag.StringVar(&jobManagerFlag, "job-manager", "", "SLURM|PBS|SGE|NONE")

	command.Register(cmd)
}

type aligners []string

// String is the method to format the flag's value, part of the flag.Value interface.
// The String method's output will be used in diagnostics.
func (a *aligners) String() string {
	return strings.Join(*a, ", ")
}

// Set is the method to set the flag value, part of the flag.Value interface.
// Set's argument is a string to be parsed to set the flag.
// It's a comma-separated list, so we split it.
func (a *aligners) Set(value string) error {
	// If we wanted to allow the flag to be set multiple times,
	// accumulating values, we would delete this if statement.
	// That would permit usages such as
	//	-deltaT 10s -deltaT 15s
	// and other combinations.
	//if len(*i) > 0 {
	//	return errors.New("interval flag already set")
	//}
	for _, aligner := range strings.Split(value, ",") {
		/*
			switch aligner {
			default:
				return ErrUnsupportedAligner{aligner: aligner}
			case "bowtie2":
				fallthrough
			case "novoalign":
				*a = append(*a, aligner)
			}
		*/
		*a = append(*a, aligner)
	}
	return nil
}

type ErrUnsupportedAligner struct {
	aligner string
}

func (e ErrUnsupportedAligner) Error() string {
	return fmt.Sprintf("nasp: unsupported aligner '%s'", e.aligner)
}

const alignTemplate = `#!/bin/bash
#PBS -t 0-{{.Tasks}}
#PBS -l nodes=1:ppn=4,mem=1gb,walltime=4:00:00
#PBS -o logs/nasp.align.out
#PBS -e logs/nasp.align.err
#set -euo pipefail

workdir=${PBS_O_WORKDIR}
cd ${workdir}

mkdir logs || true

ncpu=${PBS_NP}
index=${PBS_ARRAYID}*{{.StepSize}}
reference={{.Reference}}


{{range $index, $element := .Files}}
files[{{$index}}]={{$element}}{{end}}


# sample_name is the filename with path and extension removed
#
# Example:
#   FASTQ=/path/to/example_001.fastq.gz
#   sample_name=example_001
#
# Strip extension
sample_name=${files[$index]%%.*}
# Strip path
sample_name=${sample_name##*/}

{{define "bowtie2"}}
{{if eq .StepSize 1}}bowtie2 -p ${ncpu} --rg-id ${sample_name} --rg SM:${sample_name} -x reference -U ${files[$index]} \
{{end}}{{if eq .StepSize 2}}bowtie2 -p ${ncpu} --rg-id ${sample_name} --rg SM:${sample_name} -x reference -1 ${files[$index]} -2 ${files[$index+1]} \{{end}}
| samtools view -@ ${ncpu} -S -b -h - \
| samtools sort -@ ${ncpu} - ${sample_name}-bowtie2
samtools index ${sample_name}-bowtie2.bam
{{end}}
{{define "novoalign"}}
novoalign -d ${reference_fasta}.idx -o SAM @RG\\tID:${sample_name}\\tSM:${sample_name} -f ${files[$index]} ${files[$index+1]} \
| samtools view -@ ${ncpu} -S -b -h - \
| samtools sort -@ ${ncpu} - ${sample_name}-novo
samtools index ${sample_name}-novo.bam
{{end}}
`

type Job struct {
	Reference string
	StepSize  int
	Tasks     int
	Files     []string
}

var tmpl = template.Must(template.New("align").Parse(alignTemplate))

func runAlign(cmd *command.Command, args []string) error {
	/*
		path, err := exec.LookPath("sbatch")
		if err != nil {
			log.Println("Job")
		}
	*/

	if referenceFlag == "" {
		return errors.New("nasp.align: no reference fasta specified")
	}

	if len(alignersFlag) < 1 {
		return errors.New("nasp.align: no aligners specified")
	}

	job := Job{
		Reference: referenceFlag,
		StepSize:  2,
		Tasks:     len(args) / 2,
		Files:     args,
	}

	if err := tmpl.Execute(os.Stdout, job); err != nil {
		return fmt.Errorf("execution failed: %s", err)
	}

	for _, aligner := range alignersFlag {
		if err := tmpl.ExecuteTemplate(os.Stdout, aligner, job); err != nil {
			return err
		}
	}

	return nil
}
