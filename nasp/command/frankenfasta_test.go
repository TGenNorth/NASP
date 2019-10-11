package command

import (
	"strings"
	"testing"
)

func TestDeltaReadFrom(t *testing.T) {
	var delta delta
	if _, err := delta.ReadFrom(strings.NewReader("NUCMER")); err != ErrUnexpectedFormat {
		t.Fatalf("want:ErrUnexpectedFormat have:%s", err)
	}
}
