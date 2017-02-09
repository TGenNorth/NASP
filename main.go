package main

import (
	"bufio"
	"flag"
	"io"
	"log"
	"os"
	"strings"
	"text/template"

	"github.com/TGenNorth/NASP/command"
	_ "github.com/TGenNorth/NASP/command/export"
	_ "github.com/TGenNorth/NASP/command/matrix"
)

// TODO: document --version flag is usage template
// version is set at compile time when built with the following command:
// go build -ldflags "-X main.version=$(git rev-parse --short HEAD)"
var version string
var versionFlag bool

var commands = command.Commands

func init() {
	flag.BoolVar(&versionFlag, "version", false, "")
}

func main() {
	flag.Usage = usage
	flag.Parse()
	log.SetFlags(0)

	if versionFlag {
		log.Printf("%s", version)
		return
	}

	args := flag.Args()
	if len(args) < 1 {
		usage()
		return
	}

	if args[0] == "help" {
		help(args[1:])
		return
	}

	for _, cmd := range commands {
		if cmd.Name() == args[0] {
			cmd.Flag.Usage = func() { cmd.Usage(nil) }
			cmd.Flag.Parse(args[1:])
			args = cmd.Flag.Args()
			if err := cmd.Run(cmd, args); err != nil {
				cmd.Usage(err)
			}
			return
		}
	}

	log.Fatalf("nasp: unknown subcommand %q\nRun 'nasp help' for usage.\n", args[0])
}

var usageTemplate = `nasp is a pipeline to collect and report on statistically-relevant,
high-confidence positions in a collection of genomes with emphasis
on single nucleotide polymorphisms (SNPs).

Usage:

	nasp command [arguments]

The commands are:
{{range .}}
	{{.Name | printf "%-11s"}} {{.Short}}{{end}}

Use "nasp help [command]" for more information about a command.
`

var helpTemplate = `usage: nasp {{.UsageLine}}

{{.Long | trim}}
`

// tmpl executes the given template text on data, writing the result to w.
func tmpl(w io.Writer, text string, data interface{}) {
	t := template.Must(template.New("root").Funcs(template.FuncMap{"trim": strings.TrimSpace}).Parse(text))
	err := t.Execute(w, data)
	if err != nil {
		panic(err)
	}
}

func printUsage(w io.Writer) {
	bw := bufio.NewWriter(w)
	tmpl(bw, usageTemplate, commands)
	bw.Flush()
}

func usage() {
	printUsage(os.Stderr)
}

// help implements the 'help' command.
func help(args []string) {
	if len(args) == 0 {
		printUsage(os.Stdout)
		return
	}

	if len(args) != 1 {
		log.Fatal("usage: nasp help command\n\nToo many arguments given.")
	}

	arg := args[0]

	for _, cmd := range commands {
		if cmd.Name() == arg {
			tmpl(os.Stdout, helpTemplate, cmd)
			return
		}
	}

	log.Fatalf("Unknown help topic %#q.  Run 'nasp help'.\n", arg)
}
