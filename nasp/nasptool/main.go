package main

import (
	"bufio"
	"flag"
	"io"
	"log"
	"os"
	"text/template"

	"github.com/TGenNorth/NASP/command"
	//_ "github.com/TGenNorth/NASP/command/align"
	_ "github.com/TGenNorth/NASP/command/export"
	//_ "github.com/TGenNorth/NASP/command/index"
	_ "github.com/TGenNorth/NASP/command/matrix"
	//_ "github.com/TGenNorth/NASP/command/snpcall"
)

var commands = command.Commands

func main() {
	flag.Usage = usage
	flag.Parse()
	log.SetFlags(0)

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
		if cmd.Name() == args[0] && cmd.Runnable() {
			cmd.Flag.Usage = func() { cmd.Usage() }
			if cmd.CustomFlags {
				args = args[1:]
			} else {
				cmd.Flag.Parse(args[1:])
				args = cmd.Flag.Args()
			}
			cmd.Run(cmd, args)
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
{{range .}}{{if .Runnable}}
	{{.Name | printf "%-11s"}} {{.Short}}{{end}}{{end}}

Use "nasp help [command]" for more information about a command.

Additional help topics:
{{range .}}{{if not .Runnable}}
	{{.Name | printf "%-11s"}} {{.Short}}{{end}}{{end}}

Use "nasp help [topic]" for more information about that topic.

`

var helpTemplate = `{{if .Runnable}}usage: nasp {{.UsageLine}}

{{end}}{{.Long | trim}}
`

// tmpl executes the given template text on data, writing the result to w.
func tmpl(w io.Writer, text string, data interface{}) {
	t := template.Must(template.New("root").Parse(text))
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
