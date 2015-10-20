package command

import (
	"flag"
	"fmt"
	"os"
	"strings"
)

// Command is a nasp subcommand
type Command struct {
	// Run runs the command.
	// The args are the arguments after the command name.
	Run func(cmd *Command, args []string) error

	// UsageLine is the one-line usage message.
	// The first word in the line is taken to be the command name.
	UsageLine string

	// Short is the short description shown in the 'go help' output.
	Short string

	// Long is the long message shown in the 'go help <this-command>' output.
	Long string

	// Flag is a set of flags specific to this command.
	Flag flag.FlagSet
}

// Name returns the command's name: the first word in the usage line.
func (c *Command) Name() string {
	name := c.UsageLine
	i := strings.Index(name, " ")
	if i >= 0 {
		name = name[:i]
	}
	return name
}

func (c *Command) Usage(err error) {
	fmt.Fprintf(os.Stderr, "usage: %s\n\n", c.UsageLine)
	fmt.Fprintf(os.Stderr, "%s\n", strings.TrimSpace(c.Long))
	if err != nil {
		fmt.Fprintf(os.Stderr, "\n%s\n", err)
	}
	os.Exit(1)
}

// Commands lists the available commands and help topics.
// The order here is the order in which they are printed by 'nasp help'.
var Commands = []*Command{
	cmdDuplicates,
	cmdFrankenfasta,
}

func Register(cmd *Command) {
	Commands = append(Commands, cmd)
}
