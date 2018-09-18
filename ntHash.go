// Package ntHash is a port of ntHash (https://github.com/bcgsc/ntHash) recursive hash function for DNA kmers.
//
// It was inspired by the Rust port by Luiz Irber (https://github.com/luizirber/nthash)
//
package ntHash

import (
	"fmt"
)

// MAXIMUM_K_SIZE is the maximum k-mer size permitted by the ntHash iterator
const MAXIMUM_K_SIZE = uint(31)

// hash takes a base (byte) and returns ....
func hash(base uint8) uint64 {
	switch base {
	case 'A':
		return 0x3c8bfbb395c60474
	case 'C':
		return 0x3193c18562a02b4c
	case 'G':
		return 0x20323ed082572324
	case 'T':
		return 0x295549f54be24456
	case 'N':
		return 0
	default:
		//panic(fmt.Errorf("non a/c/t/g/n base: %c", base))
		return 0
	}
}

// rcHash takes a base (byte) and returns...
func rcHash(base uint8) uint64 {
	switch base {
	case 'A':
		return 0x295549f54be24456
	case 'C':
		return 0x20323ed082572324
	case 'G':
		return 0x3193c18562a02b4c
	case 'T':
		return 0x3c8bfbb395c60474
	case 'N':
		return 0
	default:
		//panic(fmt.Errorf("non a/c/t/g/n base: %c", base))
		return 0
	}
}

// roL is a function to bit shift to the left by "n" positions
func roL(v uint64, n uint) uint64 {
	if (n & 63) == 0 {
		return v
	}
	return (v << n) | (v >> (64 - n))
}

// roR is a function to bit shift to the right by "n" positions
func roR(v uint64, n uint) uint64 {
	if (n & 63) == 0 {
		return v
	}
	return (v >> n) | (v << (64 - n))
}

// ntf64 generates the ntHash for the forward strand of the kmer
func ntf64(seq []byte, i, k uint) uint64 {
	var hv uint64
	for i < k {
		hv = roL(hv, 1)
		hv ^= hash(seq[i])
		i++
	}
	return hv
}

// ntr64 generates the ntHash for the reverse strand of the kmer
func ntr64(seq []byte, i, k uint) uint64 {
	var hv uint64
	for i < k {
		hv = roL(hv, 1)
		hv ^= rcHash(seq[k-1-i])
		i++
	}
	return hv
}

// ntc64 generates the canonical ntHash
func ntc64(seq []byte, i, k uint) uint64 {
	fh := ntf64(seq, i, k)
	rh := ntr64(seq, i, k)
	if rh < fh {
		return rh
	}
	return fh
}

// nthash returns the canonical ntHash for each k-mer in a sequence
// it does not use the rolling hash properties of ntHash
func nthash(seq []byte, k int) []uint64 {
	hvs := make([]uint64, (len(seq) - (k - 1)))
	for i := 0; i <= (len(seq) - k); i++ {
		hvs[i] = ntc64(seq[i:i+k], 0, uint(k))
	}
	return hvs
}

// ntHash iterator struct
type nthi struct {
	seq        []byte
	k          uint
	fh         uint64
	rh         uint64
	currentIdx uint
	maxIdx     uint
}

// New creates a new ntHash iterator
func New(seq []byte, kSize int) (*nthi, error) {
	seqLen := uint(len(seq))
	k := uint(kSize)
	if k > seqLen {
		return nil, fmt.Errorf("k size is greater than sequence length (%d vs %d)!", k, seqLen)
	}
	if k > MAXIMUM_K_SIZE {
		return nil, fmt.Errorf("k size is greater than the maximum allowed k size (%d vs %d)!", k, MAXIMUM_K_SIZE)
	}
	fh := ntf64(seq[0:k], 0, k)
	rh := ntr64(seq[0:k], 0, k)
	nthi := &nthi{
		seq:        seq,
		k:          k,
		fh:         fh,
		rh:         rh,
		currentIdx: 0,
		maxIdx:     (seqLen - (k - 1)),
	}
	return nthi, nil
}

// Next method returns the next canonical ntHash value from the ntHash iterator
func (nthi *nthi) Next() uint64 {
	// end the iterator if we have got to the maximum index position TODO: this needs to be done in a better way.
	if nthi.currentIdx >= nthi.maxIdx {
		return 0
	}
	// roll the hash if index>0
	if nthi.currentIdx != 0 {
		i := nthi.currentIdx
		k := nthi.k
		// alg 3. of ntHash paper
		nthi.fh = roL(nthi.fh, 1)
		nthi.fh ^= roL(hash(nthi.seq[i-1]), nthi.k)
		nthi.fh ^= hash(nthi.seq[i+k-1])
		nthi.rh = roR(nthi.rh, 1)
		nthi.rh ^= roR(rcHash(nthi.seq[i-1]), 1)
		nthi.rh ^= roL(rcHash(nthi.seq[i+k-1]), nthi.k-1)
	}
	nthi.currentIdx++
	// return the canonical ntHash
	if nthi.rh < nthi.fh {
		return nthi.rh
	}
	return nthi.fh
}

// Hash method returns a channel to range over the canonical ntHash values of a sequence held by the ntHash iterator
func (nthi *nthi) Hash() <-chan uint64 {
	hashChan := make(chan uint64)
	go func() {
		defer close(hashChan)
		for nthi.currentIdx < nthi.maxIdx {
			hashChan <- nthi.Next()
		}
	}()
	return hashChan
}
