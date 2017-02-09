package index

import (
	"fmt"
	"os"
	"strings"
	"text/template"

	"github.com/TGenNorth/NASP/command"
)

var cmd = &command.Command{
	UsageLine: "index --aligner [aligner_name[,aligner_name,...]] --snpcaller [snpcaller[,snpcaller]] --reference [reference.fasta]",
	Short:     "create reference index files",
	Long: `
Align TODO

For more about <ACTION>, see 'nasp help <ACTION>'.
	`,
}

var (
	alignersFlag   aligners
	snpcallerFlag  snpcallers
	dryRunFlag     bool
	jobManagerFlag string
	referenceFlag  string
)

func init() {
	cmd.Run = runIndexReference
	cmd.Flag.Var(&alignersFlag, "aligner", "comma-separated list of aligners")
	cmd.Flag.Var(&snpcallerFlag, "snpcaller", "comma-separated list of snpcallers")
	cmd.Flag.BoolVar(&dryRunFlag, "dry-run", true, "write the job script, but do not run it")
	cmd.Flag.StringVar(&referenceFlag, "reference", "", "path to the reference fasta")
	cmd.Flag.StringVar(&jobManagerFlag, "job-manager", "", "SLURM|PBS|SGE|NONE")

	command.Register(cmd)
}

type aligners []string

// String is the method to format the flag's value, part of the flag.Value interface.
// The String method's output will be used in diagnostics.
func (a *aligners) String() string {
	return fmt.Sprint(*a)
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
		switch aligner {
		default:
			return ErrUnsupportedAligner{aligner: aligner}
		case "bowtie2":
			fallthrough
		case "novoalign":
			*a = append(*a, aligner)
		}
	}
	return nil
}

type snpcallers []string

// String is the method to format the flag's value, part of the flag.Value interface.
// The String method's output will be used in diagnostics.
func (s *snpcallers) String() string {
	return fmt.Sprint(*s)
}

// Set is the method to set the flag value, part of the flag.Value interface.
// Set's argument is a string to be parsed to set the flag.
// It's a comma-separated list, so we split it.
func (s *snpcallers) Set(value string) error {
	// If we wanted to allow the flag to be set multiple times,
	// accumulating values, we would delete this if statement.
	// That would permit usages such as
	//	-deltaT 10s -deltaT 15s
	// and other combinations.
	//if len(*i) > 0 {
	//	return errors.New("interval flag already set")
	//}
	for _, snpcaller := range strings.Split(value, ",") {
		switch snpcaller {
		default:
			return ErrUnsupportedSnpcaller{snpcaller: snpcaller}
		case "gatk":
			*s = append(*s, snpcaller)
		}
	}
	return nil
}

type ErrUnsupportedAligner struct {
	aligner string
}

func (e ErrUnsupportedAligner) Error() string {
	return fmt.Sprintf("nasp: unsupported aligner '%s'", e.aligner)
}

type ErrUnsupportedSnpcaller struct {
	snpcaller string
}

func (e ErrUnsupportedSnpcaller) Error() string {
	return fmt.Sprintf("nasp: unsupported snpcaller '%s'", e.snpcaller)
}

const alignTemplate = `#!/bin/bash
#PBS -t 0-{{.Tasks}}
#PBS -l nodes=1:ppn=4,mem=1gb,walltime=4:00:00
#PBS -o logs/index.reference.out
#PBS -e logs/index.reference.err
#set -euo pipefail

workdir=${PBS_O_WORKDIR}
cd ${workdir}

mkdir logs || true

ncpu=${PBS_NP}
index=${PBS_ARRAYID}*{{.StepSize}}
reference_fasta={{.Reference}}


# sample_name is the filename with path and extension removed
#
# Example:
#   files[0]=/path/to/example_001.fastq.gz
#   sample_name=example_001
#
# Strip extension
sample_name=${files[$index]%%.*}
# Strip path
sample_name=${sample_name##*/}

{{define "bowtie2"}}
bowtie2-build "${reference_fasta}" reference
{{end}}

{{define "gatk"}}
# GATK
#http://gatkforums.broadinstitute.org/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference
java -Xmx1G -jar /packages/tnorth/bin/CreateSequenceDictionary.jar R=${reference_fasta} O=${reference_fasta%.*}.dict
samtools faidx ${reference_fasta}
{{end}}

{{define "novoalign"}}
novoindex -t ${NCPU} ${REFERENCE_FASTA}.idx ${REFERENCE_FASTA}
{{end}}
`

type Job struct {
	Reference string
	StepSize  int
	Tasks     int
	Files     []string
}

var tmpl = template.Must(template.New("root").Parse(alignTemplate))

func runIndexReference(cmd *command.Command, args []string) error {
	/*
		path, err := exec.LookPath("sbatch")
		if err != nil {
			log.Println("Job")
		}
	*/

	if err := tmpl.Execute(os.Stdout, Job{
		Reference: referenceFlag,
		StepSize:  2,
		Tasks:     len(args) / 2,
		Files:     args,
	}); err != nil {
		return fmt.Errorf("execution failed: %s", err)
	}

	for _, aligner := range alignersFlag {
		if err := tmpl.ExecuteTemplate(os.Stdout, aligner, ""); err != nil {
			return err
		}
	}
	for _, snpcaller := range snpcallerFlag {
		if err := tmpl.ExecuteTemplate(os.Stdout, snpcaller, ""); err != nil {
			return err
		}
	}

	return nil

}
