package command

import (
	"strings"
	"testing"
)

func TestDeltaReadFrom(t *testing.T) {
	d := &delta{}

	// encountered when the query assembly was an empty file
	if _, err := d.ReadFrom(strings.NewReader("NUCMER")); err != ErrUnexpectedFormat {
		t.Fatalf("want:ErrUnexpectedFormat have:%s", err)
	}

	d = &delta{}
	// encountered from a bbtools generated dataset:
	//   randomgenome.sh len=1000 out=reference.fasta
	//   mutate.sh id=0.9 in=reference.fasta out=query.fasta
	r := strings.NewReader("/path/to/reference.fasta /path/to/query.fasta\nNUCMER")
	if _, err := d.ReadFrom(r); err != ErrUnexpectedFormat {
		t.Fatalf("want:ErrUnexpectedFormat have:%v", err)
	}
}
