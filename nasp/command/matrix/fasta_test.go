package matrix

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"testing"
)

func tmpFasta(t *testing.T) *os.File {
	file, err := ioutil.TempFile("", "")
	if err != nil {
		t.Fatal(err)
	}
	if _, err := file.Write([]byte(">ContigA\nGATC\nABCD\r\nefgh\n>ContigB\nCTAGDCBA")); err != nil {
		t.Fatal(err)
	}
	return file
}

func closeFasta(file *os.File) {
	file.Close()
	os.Remove(file.Name())
}

func TestNewReferenceInvalidReference(t *testing.T) {
	refFasta := tmpFasta(t)
	defer closeFasta(refFasta)

	_, err := NewReference("invalid", "")
	if err == nil {
		t.Fatal(err)
	}

	t.Log(err.Error())
}

func TestReferenceDups(t *testing.T) {
	if testing.Short() {
		t.Skip()
	}

	refFilepath := "./testdata/FmcMrsaReference.txt"
	dupFilepath := "./testdata/FmcMrsaDuplicates.txt"

	file, err := os.Open(dupFilepath)
	if err != nil {
		t.Fatal(err)
	}
	br := bufio.NewReader(file)

	reference, err := NewReference(refFilepath, dupFilepath)
	if err != nil {
		t.Fatal(err)
	}

	var dup []byte
	for {
		name, err := reference.NextContig()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}

		line, _, err := br.ReadLine()
		if err != nil {
			t.Fatal(err)
		}

		if !bytes.Equal([]byte(name), line[1:]) {
			t.Errorf("Contig name should match - expect: %s actual: `%s`\n", line[1:], name)
			t.Log(reference.ref.isPrefix, reference.dup.isPrefix, len(name))
		}

		for isPrefix := true; isPrefix; {
			// TODO: try with a static 80
			// len(line) should be a static 80 characters.
			_, dup, isPrefix, err = reference.ReadPositions(80)
			if err != nil {
				t.Fatal(err)
			}

			if !isPrefix {
				// TODO: Check dup is empty
				break
			}

			line, _, err := br.ReadLine()
			if err != nil {
				t.Fatal(err)
			}

			if !bytes.Equal(line, dup) {
				t.Errorf("Expect:\n%s\nActual:\n%s\n", line, dup)
				t.Log("dup", dup, "len(dup)", len(dup), reference.dup)
				var offset int64
				for {
					line, err := br.ReadSlice('>')
					offset += int64(len(line))
					if err == nil {
						break
					}
					if err != bufio.ErrBufferFull {
						t.Fatal(err)
					}
				}
				t.Logf("Found a contig marker (>) after reading %d characters\n", offset)
				/*
					line, isPrefix, err := br.ReadLine()
					if err != nil {
						t.Fatal(err)
					}
					if isPrefix {
						t.Fatalf("Failed to read contig description into buffer: %s\n", line)
					}
					t.Logf("Found contig: %s\n", line)
				*/
				br.UnreadByte()
			}
		}
	}

}

func TestCanRebuildReferenceAndDupsFromParserOutput(t *testing.T) {
	if testing.Short() {
		t.Skip()
	}

	// Given a fasta file formatted to 80 char lines,
	// we should be able to reconstruct the file from the Fasta parser calls.
	refFilepath := "./testdata/FmcMrsaReference.txt"
	dupFilepath := "./testdata/FmcMrsaDuplicates.txt"

	expectRef, err := ioutil.ReadFile(refFilepath)
	if err != nil {
		t.Fatal(err)
	}

	expectDup, err := ioutil.ReadFile(dupFilepath)
	if err != nil {
		t.Fatal(err)
	}

	actualRef := bytes.NewBuffer(make([]byte, 0, len(expectRef)))
	actualDup := bytes.NewBuffer(make([]byte, 0, len(expectDup)))

	reference, err := NewReference(refFilepath, dupFilepath)
	if err != nil {
		t.Fatal(err)
	}

	for {
		name, err := reference.NextContig()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}
		actualRef.WriteByte('>')
		actualDup.WriteByte('>')
		actualRef.WriteString(name)
		actualDup.WriteString(name)
		actualRef.WriteByte('\n')
		actualDup.WriteByte('\n')

		var ref, dup []byte
		for isPrefix := true; isPrefix; {
			ref, dup, isPrefix, err = reference.ReadPositions(80)
			if err != nil {
				t.Error(err)
			}
			actualRef.Write(ref)
			actualDup.Write(dup)
			actualRef.WriteByte('\n')
			actualDup.WriteByte('\n')
		}
	}

	if !bytes.Equal(expectRef, actualRef.Bytes()) {
		t.Error("The parsed reference does not match the source reference")
	}
	if !bytes.Equal(expectDup, actualDup.Bytes()) {
		t.Error("The parsed dups does not match the source dups")
	}
}

func TestDup(t *testing.T) {
	if testing.Short() {
		t.Skip()
	}

	// Given a fasta file formatted to 80 char lines,
	// we should be able to reconstruct the file from the Fasta parser calls.
	dupsFilepath := "./testdata/FmcMrsaDuplicates.txt"

	expected, err := ioutil.ReadFile(dupsFilepath)
	if err != nil {
		t.Fatal(err)
	}

	actual := bytes.NewBuffer(make([]byte, 0, len(expected)))

	dup, err := NewFasta(&NaspFile{Filepath: dupsFilepath}, false)
	if err != nil {
		t.Fatal(err)
	}

	for {
		name, err := dup.NextContig()
		if err == io.EOF {
			break
		}
		if err != nil {
			t.Fatal(err)
		}

		actual.WriteByte('>')
		actual.WriteString(name)
		actual.WriteByte('\n')

		for dup.isPrefix {
			calls, _, _, err := dup.ReadPositions(80)
			if err != nil {
				t.Fatal(err)
			}
			actual.Write(calls)
			actual.WriteByte('\n')
		}
	}

	if !bytes.Equal(expected, actual.Bytes()) {
		t.Fatal("Result does not match")
	}
}

func TestReference(t *testing.T) {
	t.Skip("This test prints out the rebuilt contents of a duplicates file from the Reference parser output")
	ref, err := NewReference("./test_time/reference/reference.fasta", "./test_time/reference/duplicates.txt")
	if err != nil {
		t.Fatal(err)
	}
	for {
		name, _ := ref.NextContig()
		fmt.Println(">", name)

		for isPrefix := true; isPrefix; {
			var dup []byte
			_, dup, isPrefix, err = ref.ReadPositions(80)
			if err != nil {
				t.Fatal(err)
			}
			fmt.Printf("%s", dup)
		}
		fmt.Println()
	}
}

func TestNewReferenceInvalidDuplicates(t *testing.T) {
	refFasta := tmpFasta(t)
	defer closeFasta(refFasta)

	dupFasta := tmpFasta(t)
	defer closeFasta(dupFasta)

	_, err := NewReference(refFasta.Name(), "invalid")
	if err == nil {
		t.Fatal(err)
	}

	t.Log(err.Error())
}

func TestNewReferenceOnly(t *testing.T) {
	refFasta := tmpFasta(t)
	defer closeFasta(refFasta)

	if ref, err := NewReference(refFasta.Name(), ""); ref == nil || err != nil {
		t.Fatal(err)
	}
}

func TestNewReferenceWithDuplicates(t *testing.T) {
	refFasta := tmpFasta(t)
	defer closeFasta(refFasta)

	dupFasta := tmpFasta(t)
	defer closeFasta(dupFasta)

	if ref, err := NewReference(refFasta.Name(), dupFasta.Name()); ref == nil || err != nil {
		t.Fatal(err)
	}
}

func TestNext(t *testing.T) {
	file := tmpFasta(t)
	defer closeFasta(file)

	fasta, err := NewFasta(&NaspFile{Filepath: file.Name()}, false)
	if err != nil {
		t.Fatal(err)
	}

	expect := "ContigA"
	name, err := fasta.NextContig()
	if name != expect || err != nil {
		t.Fatalf("Expected: %s, nil Observed: %s, %s\n", name, err)
	}

	expect = "ContigB"
	name, err = fasta.NextContig()
	if name != expect || err != nil {
		t.Fatalf("Expected: %s, nil Observed: %s, %s\n", name, err)
	}

	expect = ""
	name, err = fasta.NextContig()
	if name != expect || err != io.EOF {
		t.Fatalf("Expected: %s, io.EOF Observed: %s, %s\n", name, err)
	}
}

/*
func TestIndex(t *testing.T) {
	file, err := ioutil.TempFile("", "")
	if err != nil {
		t.Fatal(err)
	}
	defer func() {
		file.Close()
		os.Remove(file.Name())
	}()

	if _, err := file.Write([]byte(">ContigA\nGATC\nABCD\r\nefgh\r\n>ContigB\nCTAGDCBA")); err != nil {
		t.Fatal(err)
	}

	fasta, err := NewFasta(file.Name(), true)
	if err != nil {
		t.Fatal(err)
	}

	done := make(chan struct{})
	defer close(done)

	contigName := "ContigA"
	if err = fasta.SeekContig(contigName); err != nil {
		t.Fatal(err)
	}
	expect := []byte("GATC")
	actual, err := fasta.ReadPositions(4)
	if !bytes.Equal(expect, actual) {
		t.Fatalf("%s expected: %s actual: %s", contigName, expect, actual)
	}
	expect = []byte("ABCDefgh")
	actual, err = fasta.ReadPositions(8)
	if !bytes.Equal(expect, actual) {
		t.Fatalf("%s expected: %s actual: %s", contigName, expect, actual)
	}
	// TODO: Expected it yields empty positions when the contig is exhausted

	contigName = "ContigB"
	if err = fasta.SeekContig(contigName); err != nil {
		t.Fatal(err)
	}
	expect = []byte("CTAGDCBA")
	actual, err = fasta.ReadPositions(8)
	if !bytes.Equal(expect, actual) {
		t.Fatalf("%s expected: %c actual: %c", contigName, expect, actual)
	}

	// TODO: Expected it yields empty positions when the file is exhausted
	// until the channel is closed
}
*/
