// The subcommand design pattern is based on the Go command.
// All other source code is copyright TGenNorth.

// Copyright 2011 The Go Authors.  All rights reserved.
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.
package main

import (
	"bufio"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"os/exec"
	"strings"
	"sync"
	"text/template"
	"unicode"
	"unicode/utf8"

	"github.com/TGenNorth/nasp/command"
	//_ "github.com/TGenNorth/nasp/command/align"
	_ "github.com/TGenNorth/nasp/command/export"
	//_ "github.com/TGenNorth/nasp/command/index"
	_ "github.com/TGenNorth/nasp/command/matrix"
	//_ "github.com/TGenNorth/nasp/command/snpcall"
)

var commands = command.Commands
var exitStatus = 0
var exitMu sync.Mutex

func setExitStatus(n int) {
	exitMu.Lock()
	if exitStatus < n {
		exitStatus = n
	}
	exitMu.Unlock()
}

func main() {
	flag.Usage = usage
	flag.Parse()
	log.SetFlags(0)

	args := flag.Args()
	if len(args) < 1 {
		usage()
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
			exit()
			return
		}
	}

	fmt.Fprintf(os.Stderr, "nasp: unknown subcommand %q\nRun 'nasp help' for usage.\n", args[0])
	setExitStatus(2)
	exit()
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

// An errWriter wraps a writer, recording whether a write error occurred.
type errWriter struct {
	w   io.Writer
	err error
}

func (w *errWriter) Write(b []byte) (int, error) {
	n, err := w.w.Write(b)
	if err != nil {
		w.err = err
	}
	return n, err
}

// tmpl executes the given template text on data, writing the result to w.
func tmpl(w io.Writer, text string, data interface{}) {
	t := template.New("top")
	t.Funcs(template.FuncMap{"trim": strings.TrimSpace, "capitalize": capitalize})
	template.Must(t.Parse(text))
	ew := &errWriter{w: w}
	err := t.Execute(ew, data)
	if ew.err != nil {
		// I/O error writing. Ignore write on closed pipe.
		if strings.Contains(ew.err.Error(), "pipe") {
			os.Exit(1)
		}
		fatalf("writing output: %v", ew.err)
	}
	if err != nil {
		panic(err)
	}
}

func capitalize(s string) string {
	if s == "" {
		return s
	}
	r, n := utf8.DecodeRuneInString(s)
	return string(unicode.ToTitle(r)) + s[n:]
}

func printUsage(w io.Writer) {
	bw := bufio.NewWriter(w)
	tmpl(bw, usageTemplate, commands)
	bw.Flush()
}

func usage() {
	printUsage(os.Stderr)
	os.Exit(2)
}

// help implements the 'help' command.
func help(args []string) {
	if len(args) == 0 {
		printUsage(os.Stdout)
		// not exit 2: succeeded at 'nasp help'.
		return
	}
	if len(args) != 1 {
		fmt.Fprintf(os.Stderr, "usage: nasp help command\n\nToo many arguments given.\n")
		os.Exit(2) // failed at 'nasp help'
	}

	arg := args[0]

	for _, cmd := range commands {
		if cmd.Name() == arg {
			tmpl(os.Stdout, helpTemplate, cmd)
			// not exit 2: succeeded at 'nasp help cmd'.
			return
		}
	}

	fmt.Fprintf(os.Stderr, "Unknown help topic %#q.  Run 'nasp help'.\n", arg)
	os.Exit(2) // failed at 'nasp help cmd'
}

var atexitFuncs []func()

func atexit(f func()) {
	atexitFuncs = append(atexitFuncs, f)
}

func exit() {
	for _, f := range atexitFuncs {
		f()
	}
	os.Exit(exitStatus)
}

func fatalf(format string, args ...interface{}) {
	errorf(format, args...)
	exit()
}

func errorf(format string, args ...interface{}) {
	log.Printf(format, args...)
	setExitStatus(1)
}

var logf = log.Printf

func exitIfErrors() {
	if exitStatus != 0 {
		exit()
	}
}

func run(cmdargs ...interface{}) {
	cmdline := stringList(cmdargs...)
	cmd := exec.Command(cmdline[0], cmdline[1:]...)
	cmd.Stdout = os.Stdout
	cmd.Stderr = os.Stderr
	if err := cmd.Run(); err != nil {
		errorf("%v", err)
	}
}

func runOut(dir string, cmdargs ...interface{}) []byte {
	cmdline := stringList(cmdargs...)
	cmd := exec.Command(cmdline[0], cmdline[1:]...)
	cmd.Dir = dir
	out, err := cmd.CombinedOutput()
	if err != nil {
		os.Stderr.Write(out)
		errorf("%v", err)
		out = nil
	}
	return out
}

// stringList's arguments should be a sequence of string or []string values.
// stringList flattens them into a single []string.
func stringList(args ...interface{}) []string {
	var x []string
	for _, arg := range args {
		switch arg := arg.(type) {
		case []string:
			x = append(x, arg...)
		case string:
			x = append(x, arg)
		default:
			panic("stringList: invalid argument of type " + fmt.Sprintf("%T", arg))
		}
	}
	return x
}
