package snpcall

import (
	"errors"
	"fmt"
	"log"
	"os"
	"strings"
	"text/template"

	"github.com/TGenNorth/nasp/command"
)

var cmd = &command.Command{
	UsageLine: "snpcall --snpcaller [snpcaller_name[,snpcaller_name,...]] --reference [reference.fasta]",
	Short:     "aggregate sample calls into matrices and stats",
	Long: `
Align TODO

For more about <ACTION>, see 'nasp help <ACTION>'.
	`,
}

var (
	snpcallersFlag snpcallers
	dryRunFlag     bool
	jobManagerFlag string
	referenceFlag  string
)

func init() {
	cmd.Run = runSnpcall
	cmd.Flag.Var(&snpcallersFlag, "snpcaller", "comma-separated list of snpcallers")
	cmd.Flag.BoolVar(&dryRunFlag, "dry-run", true, "write the job script, but do not run it")
	cmd.Flag.StringVar(&referenceFlag, "reference", "", "path to the reference fasta")
	cmd.Flag.StringVar(&jobManagerFlag, "job-manager", "", "SLURM|PBS|SGE|NONE")

	command.Register(cmd)
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
		case "varscan":
			fallthrough
		case "gatk":
			*s = append(*s, snpcaller)
		}
	}
	return nil
}

type ErrUnsupportedSnpcaller struct {
	snpcaller string
}

func (e ErrUnsupportedSnpcaller) Error() string {
	return fmt.Sprintf("nasp: unsupported aligner '%s'", e.snpcaller)
}

const snpcallTemplate = `#!/bin/bash
#PBS -t 0-{{.Tasks}}
#PBS -l nodes=1:ppn=4,mem=1gb,walltime=4:00:00
#PBS -o logs/bowtie2.align.out
#PBS -e logs/bowtie2.align.err
#set -euo pipefail

workdir=${PBS_O_WORKDIR}
cd ${workdir}

mkdir logs || true

ncpu=${PBS_NP}
index=${PBS_ARRAYID}*{{.StepSize}}
reference_fasta={{.Reference}}


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

bam=${files[$index]}

{{define "gatk"}}
java -Xmx1G -jar /packages/tnorth/bin/GenomeAnalysisTK.jar -T UnifiedGenotyper -dt NONE -glm BOTH -I ${bam} -R ${reference_fasta} -nt ${ncpu} -o ${sample_name}-gatk.vcf -out_mode EMIT_ALL_CONFIDENT_SITES -baq RECALCULATE
{{end}}
{{define "varscan"}}
samtools mpileup -B -d 10000000 -f ${reference_fasta} ${bam} |
java -Xmx1G -jar /packages/tnorth/bin/VarScan.jar mpileup2cns - --output-vcf 1 --vcf-sample-list <(echo ${sample_name}) > ${sample_name}-varscan.vcf
{{end}}
`

type Job struct {
	Reference string
	StepSize  int
	Tasks     int
	Files     []string
}

var tmpl = template.Must(template.New("root").Parse(snpcallTemplate))

func runSnpcall(cmd *command.Command, args []string) error {
	/*
		path, err := exec.LookPath("sbatch")
		if err != nil {
			log.Println("Job")
		}
	*/

	if len(snpcallersFlag) < 1 {
		log.Println("nasp.snpcall: no snpcallers listed")
	}

	if len(args) < 1 {
		return errors.New("nasp.snpcall: expected a list or glob pattern of read files")
	}

	if err := tmpl.Execute(os.Stdout, Job{
		Reference: referenceFlag,
		StepSize:  1,
		Tasks:     len(args) - 1,
		Files:     args,
	}); err != nil {
		return fmt.Errorf("execution failed: %s", err)
	}

	for _, snpcaller := range snpcallersFlag {
		if err := tmpl.ExecuteTemplate(os.Stdout, snpcaller, ""); err != nil {
			return err
		}
	}

	return nil
}
