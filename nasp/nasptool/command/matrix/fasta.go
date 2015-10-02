package matrix

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
)

type Fasta struct {
	//path string
	// TODO: aligner string
	// index maps contigName -> filePosition of the first position in each contig.
	name     string
	index    map[string]int64
	rd       *os.File
	br       *bufio.Reader
	isPrefix bool
	buf      *bytes.Buffer
}

func NewFasta(naspFile NaspFile, indexContigs bool) (*Fasta, error) {
	file, err := os.Open(naspFile.Filepath)
	if err != nil {
		return nil, err
	}

	fasta := &Fasta{
		//path:    path,
		name: naspFile.Name,
		rd:   file,
		br:   bufio.NewReader(file),
		buf:  bytes.NewBuffer(make([]byte, 0, defaultBufSize)),
	}

	if indexContigs {
		// TODO: measure memory/runtime of map vs sorted list.
		fasta.index = make(map[string]int64)
		if err = fasta.indexContigs(); err != nil {
			return nil, err
		}
	}

	return fasta, nil
}

func (f Fasta) Name() string {
	return f.name
}

/**
 * NextContig advances the io.Reader to the first position of the next contig
 * and returns its name.
 *
 * TODO: document name excludes description and prefix
 * >frankenfasta::name description
 * gatc...
 */
func (f *Fasta) NextContig() (name string, err error) {
	f.buf.Reset()
	f.isPrefix = true
	for {
		_, err = f.br.ReadSlice('>')
		switch err {
		default:
			return "", err
		case bufio.ErrBufferFull:
			break
		case nil:
			line, isPrefix, err := f.br.ReadLine()
			if err != nil {
				return "", err
			}
			if isPrefix {
				// FIXME: The default buffer size is 4096 bytes. We are not handling this situation.
				return "", fmt.Errorf("The contig description is too long to fit in the buffer: %s...\n", line)
			}
			if idx := bytes.IndexByte(line, ' '); idx != -1 {
				// Discard the description following the contig name.
				return string(line[:idx]), nil
			}
			return string(line), nil
		}
	}
}

// TODO: Indexing detect consecutive contig descriptions

// ReadPositions returns the next `n` positions or nil when the contig is
// exhausted. `n` can be any arbitrary buffer size.
func (f *Fasta) ReadPositions(n int) ([]byte, []int, []float64, error) {
	var line, peek []byte
	var err error
	if f.isPrefix {
		for f.buf.Len() < n {
			if peek, err = f.br.Peek(1); err != nil || peek[0] == '>' {
				f.isPrefix = false
				break
			}
			if line, _, err = f.br.ReadLine(); err != nil {
				line = nil
				break
			}
			// err is always nil
			f.buf.Write(line)
		}
	}
	if err == io.EOF && f.buf.Len() > 0 {
		// The file is empty, but the buffer is not.
		err = nil
	}
	return f.buf.Next(n), nil, nil, err
}

// SeekContig moves the
func (f *Fasta) SeekContig(name string) (err error) {
	if filePosition, ok := f.index[name]; ok {
		// ReadPositions will yield values from the new contig.
		f.isPrefix = true
		_, err = f.rd.Seek(filePosition, os.SEEK_SET)
		// Purge any old data from the previous file position
		// TODO: The contigs will probably be typically sequential, in which
		// case it might be efficient to advance the buffers instead of purging
		// them.
		f.br.Reset(f.rd)
		f.buf.Reset()
	} else {
		f.isPrefix = false
	}
	return err
}

/**
 * indexContigs maps the starting file position of each contig in the file.
 */
func (f Fasta) indexContigs() error {
	// TODO: Include length
	var line []byte
	var err error
	var position int

	for err != io.EOF {
		line, err = f.br.ReadSlice('>')
		switch err {
		default:
			return err
		case bufio.ErrBufferFull, io.EOF:
			position += len(line)
		case nil:
			position += len(line)
			// NOTE: ing ReadSlice could save a few megabytes of RAM over
			//	ReadBytes, but could fail with a ErrBufferFull
			if line, err = f.br.ReadBytes('\n'); err != nil {
				log.Fatal(err)
			}
			position += len(line)

			line = bytes.TrimPrefix(line, []byte("franken::"))

			var name string
			if idx := bytes.IndexAny(line, " \n"); idx < 0 {
				name = string(line)
			} else {
				name = string(line[:idx])
			}
			f.index[name] = int64(position)
		}
	}

	return nil
}

type Reference struct {
	lastContig string
	ref        *Fasta
	dup        *Fasta
}

func NewReference(refPath, dupPath string) (*Reference, error) {
	var ref, dup *Fasta
	var err error

	ref, err = NewFasta(NaspFile{Name: "Reference", Filepath: refPath}, false)
	if err != nil {
		return nil, fmt.Errorf("%s: %s", err.Error(), refPath)
	}

	if dupPath != "" {
		// duplicates.txt is optional
		// NOTE: The Name field is useful only for error messages, if at all.
		dup, err = NewFasta(NaspFile{Name: "FIXME: duplicates.txt errors should use the filename, not the Name field", Filepath: dupPath}, false)
		if err != nil {
			return nil, fmt.Errorf("%s: %s", err.Error(), dupPath)
		}
	}

	return &Reference{
		ref: ref,
		dup: dup,
	}, nil
}

// NextContig moves the Reader to the first position of the next contig and
// returns its name.
func (r *Reference) NextContig() (string, error) {
	var rname, dname string
	var rerr, derr error

	rname, rerr = r.ref.NextContig()

	// The duplicates file is optional. If absent, all positions will be
	// assumed unique.
	if r.dup == nil {
		return rname, rerr
	}

	// Typically there should be a corresponding contig for every contig in the
	// reference. If the contig is short enough, nucmer may not report anything.
	if r.lastContig == "" {
		if dname, derr = r.dup.NextContig(); derr != nil {
			log.Printf("TODO: duplicates file error handling - %s\n", derr.Error())
		}
		if dname != rname {
			// reference.fasta: ContigA duplicates.txt: ContigB
			log.Printf("The duplicates contig does not correspond to the reference contig.\n%s contig: %s %s contig: %s\n", r.ref.Name(), rname, r.dup.Name(), dname)
			r.lastContig = dname
			r.dup.isPrefix = false
		}
	} else if rname == r.lastContig {
		r.lastContig = ""
		r.dup.isPrefix = true
	} else {
		log.Printf("%s: %s %s: %s\n", r.ref.Name(), rname, r.dup.Name(), dname)
	}

	/*
		if derr != nil && derr != rerr {
			return "", fmt.Errorf("duplicates file: %s reference file: %s\n", derr.Error(), rerr.Error())
		}

		if rname != dname {
			return "", fmt.Errorf("The duplicates file should have a corresponding contig for every contig in the reference.\nThe following contigs were found in corresponding positions of the reference and duplicates files: `%s` and `%s`\n", rname, dname)
		}
	*/

	return rname, rerr
}

func (r Reference) ReadPositions(n int) (ref, dup []byte, isPrefix bool, rerr error) {
	var derr error

	ref, _, _, rerr = r.ref.ReadPositions(n)

	if r.dup != nil {
		dup, _, _, derr = r.dup.ReadPositions(n)

		if derr != nil && derr != rerr {
			return nil, nil, false, fmt.Errorf("duplicates file: %s reference file: %s\n", derr.Error(), rerr.Error())
		}
	}

	// Copy the ref and dup so they may be passed to a goroutine
	// without triggering a race condition.
	//buf := make([]byte, len(ref)+len(dup))
	buf := make([]byte, 2*len(ref))
	copy(buf[:len(ref)], ref)
	copy(buf[len(ref):], dup)
	return buf[:len(ref)], buf[len(ref):], r.ref.isPrefix, rerr
}
