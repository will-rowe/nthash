// test values have been lifted from Luiz Irber -- all credit and my thanks to him!
// see https://github.com/luizirber/nthash/blob/master/src/lib.rs
package nthash

import (
	"fmt"
	"testing"
)

var (
	kmer     = []byte("TGCAG")
	sequence = []byte("ACGTCGTCAGTCGATGCAGT")
	kmer2    = []byte("ACTGC")
)

// test seed lookup
func TestSeedLookup(t *testing.T) {
	if seedTab[kmer[0]] != 0x295549f54be24456 {
		t.Fatal()
	}
	if seedTab[kmer[1]] != 0x20323ed082572324 {
		t.Fatal()
	}
	if seedTab[kmer[2]] != 0x3193c18562a02b4c {
		t.Fatal()
	}
	if seedTab[kmer[3]] != 0x3c8bfbb395c60474 {
		t.Fatal()
	}
}

// test forward ntHash
func TestNTF64hash(t *testing.T) {
	hv := ntf64(kmer, 0, 5)
	t.Log(fmt.Printf("%x\n", hv))
	if hv != 0xbafa6728fc6dabf {
		t.Fatal()
	}
}

// test reverse ntHash
func TestNTR64(t *testing.T) {
	hv := ntr64(kmer, 0, 5)
	t.Log(fmt.Printf("%x\n", hv))
	if hv != 0x8cf2d4072cca480e {
		t.Fatal()
	}
}

// test the canonical ntHash
func TestNTC64(t *testing.T) {
	hv := ntc64(kmer, 0, 5)
	t.Log(fmt.Printf("%x\n", hv))
	if hv != 0xbafa6728fc6dabf {
		t.Fatal()
	}
}

// test the ntHash function TODO: actually test this....
func TestNTHash(t *testing.T) {
	hvs := nthash(sequence, 5)
	for i, h := range hvs {
		t.Log(i, h)
	}
}

// test the ntHash iterator constructor
func TestNewHasherNTHI(t *testing.T) {
	if _, err := NewHasher(&kmer, 10); err == nil {
		t.Fatal("should trigger k > seq error")
	}
	if _, err := NewHasher(&kmer, 200); err == nil {
		t.Fatal("should trigger k > max_k error")
	}
	nthi, err := NewHasher(&sequence, 5)
	if err != nil {
		t.Fatal()
	}
	t.Log(nthi)
}

// test the ntHash iterator next method
func TestNext(t *testing.T) {
	nthi, err := NewHasher(&kmer2, 3)
	if err != nil {
		t.Fatal()
	}
	// should return the pre-calculated ntHash for the first canonical k-mer

	if h, _ := nthi.Next(true); h != 0x9b1eda9a185413ce {
		t.Fatal()
	}
	t.Log(nthi)
	// should calculate the next canonical k-mer ntHash and return it
	if h, _ := nthi.Next(true); h != 0x9f6acfa2235b86fc {
		t.Fatal()
	}
	// should calculate the final canonical k-mer ntHash and return it
	if h, _ := nthi.Next(true); h != 0xd4a29bf149877c5c {
		t.Fatal()
	}
}

// test the ntHash iterator hash method
func TestHash(t *testing.T) {
	nthi, err := NewHasher(&kmer2, 3)
	if err != nil {
		t.Fatal()
	}
	counter := 0
	// use the canonical switch
	for hash := range nthi.Hash(true) {
		t.Log(hash)
		counter++
		switch counter {
		case 1:
			if hash != 0x9b1eda9a185413ce {
				t.Fatal()
			}
		case 2:
			if hash != 0x9f6acfa2235b86fc {
				t.Fatal()
			}
		case 3:
			if hash != 0xd4a29bf149877c5c {
				t.Fatal()
			}
		default:
			t.Fatal("unexpected output from nthi")
		}
	}
	if counter != 3 {
		t.Fatal("wrong iteration")
	}
}

// test the ntHash iterator multihash method
func TestMultiHash(t *testing.T) {
	nthi, err := NewHasher(&kmer2, 3)
	if err != nil {
		t.Fatal()
	}
	counter := 0

	// use the canonical switch and 3 multihashes
	for hashes := range nthi.MultiHash(true, 3) {
		t.Log(hashes)
		counter++
		switch counter {
		case 1:
			if hashes[0] != 0x9b1eda9a185413ce {
				t.Fatal()
			}
		case 2:
			if hashes[0] != 0x9f6acfa2235b86fc {
				t.Fatal()
			}
		case 3:
			if hashes[0] != 0xd4a29bf149877c5c {
				t.Fatal()
			}
		default:
			t.Fatal("unexpected output from nthi")
		}
	}
	if counter != 3 {
		t.Fatal("wrong iteration")
	}
}

// run benchmarks of ntHash
func BenchmarkHash(b *testing.B) {
	// run the ntHash iterator b.N times
	for n := 0; n < b.N; n++ {
		nthi, err := NewHasher(&sequence, 7)
		if err != nil {
			b.Fatal()
		}
		for range nthi.Hash(false) {
		}
	}
}

func BenchmarkCanonicalHash(b *testing.B) {
	// run the ntHash iterator b.N times
	for n := 0; n < b.N; n++ {
		nthi, err := NewHasher(&sequence, 7)
		if err != nil {
			b.Fatal()
		}
		for range nthi.Hash(true) {
		}
	}
}
