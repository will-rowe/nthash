// Package nthash is a port of ntHash (https://github.com/bcgsc/ntHash) recursive hash function for DNA kmers.
//
// It was inspired by the Rust port by Luiz Irber (https://github.com/luizirber/nthash)
//
package nthash

import (
	"fmt"
	"math"
	"sync"
)

const (
	// maxK is the maximum k-mer size permitted
	maxK uint = math.MaxUint32

	// bufferSize is the size of te buffer used by the channel in the Hash method
	bufferSize uint = 128

	// offset is used as a mask to retrieve a base's complement in the seed table
	offset uint8 = 0x07

	// seedA is the 64-bit random seed corresponding to base A
	seedA uint64 = 0x3c8bfbb395c60474

	// seedC is the 64-bit random seed corresponding to base C
	seedC uint64 = 0x3193c18562a02b4c

	// seedG is the 64-bit random seed corresponding to base G
	seedG uint64 = 0x20323ed082572324

	// seedT is the 64-bit random seed corresponding to base T
	seedT uint64 = 0x295549f54be24456

	// seedN is the 64-bit random seed corresponding to N
	seedN uint64 = 0x0000000000000000

	// seed for gerenerating multiple hash values
	multiSeed uint64 = 0x90b45d39fb6da1fa

	// multiShift is used for gerenerating multiple hash values
	multiShift uint = 27
)

// seedTab is the lookup table for the bases on their complements
var seedTab = [256]uint64{
	seedN, seedT, seedN, seedG, seedA, seedA, seedN, seedC, // 0..7
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 8..15
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 16..23
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 24..31
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 32..39
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 40..47
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 48..55
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 56..63
	seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 64..71
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 72..79
	seedN, seedN, seedN, seedN, seedT, seedT, seedN, seedN, // 80..87
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 88..95
	seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 96..103
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 104..111
	seedN, seedN, seedN, seedN, seedT, seedT, seedN, seedN, // 112..119
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 120..127
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 128..135
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 136..143
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 144..151
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 152..159
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 160..167
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 168..175
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 176..183
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 184..191
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 192..199
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 200..207
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 208..215
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 216..223
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 224..231
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 232..239
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 240..247
	seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 248..255
}

// NTHi is the ntHash iterator
type NTHi struct {
	seq        *[]byte // the sequence being hashed
	k          uint    // the k-mer size
	fh         uint64  // the current forward hash value
	rh         uint64  // the current reverse hash value
	currentIdx uint    // the current index position in the sequence being hashed
	maxIdx     uint    // the maximum index position to hash up to
}

// use object pool to reducing GC load for computation of huge number of sequences.
var poolNTHi = &sync.Pool{New: func() interface{} {
	return &NTHi{}
}}

// NewHasher is the constructor function for the ntHash iterator
// seq is a pointer to the sequence being hashed
// k is the k-mer size to use
func NewHasher(seq *[]byte, k uint) (*NTHi, error) {
	seqLen := uint(len(*seq))
	if k > seqLen {
		return nil, fmt.Errorf("k size is greater than sequence length (%d vs %d)", k, seqLen)
	}
	if k > maxK {
		return nil, fmt.Errorf("k size is greater than the maximum allowed k size (%d vs %d)", k, maxK)
	}
	fh := ntf64((*seq)[0:k], 0, k)
	rh := ntr64((*seq)[0:k], 0, k)

	nthi := poolNTHi.Get().(*NTHi)
	nthi.seq = seq
	nthi.k = k
	nthi.fh = fh
	nthi.rh = rh
	nthi.currentIdx = 0
	nthi.maxIdx = seqLen - (k - 1)

	return nthi, nil
}

// Next returns the next ntHash value from an ntHash iterator
func (nthi *NTHi) Next(canonical bool) (uint64, bool) {

	// end the iterator if we have got to the maximum index position TODO: this needs to be done in a better way.
	if nthi.currentIdx >= nthi.maxIdx {
		poolNTHi.Put(nthi)
		return 0, false
	}

	// roll the hash if index>0
	if nthi.currentIdx != 0 {
		prevBase := (*nthi.seq)[nthi.currentIdx-1]
		endBase := (*nthi.seq)[nthi.currentIdx+nthi.k-1]
		// alg 3. of ntHash paper
		nthi.fh = roL(nthi.fh, 1)
		nthi.fh ^= roL(seedTab[prevBase], nthi.k)
		nthi.fh ^= seedTab[endBase]
		nthi.rh = roR(nthi.rh, 1)
		nthi.rh ^= roR(seedTab[prevBase&offset], 1)
		nthi.rh ^= roL(seedTab[endBase&offset], nthi.k-1)
	}
	nthi.currentIdx++

	if canonical {
		return nthi.getCanonical(), true
	}
	return nthi.fh, true
}

// Hash returns a channel to range over the canonical ntHash values of a sequence
// canonical is set true to return the canonical k-mers, otherwise the forward hashes are returned
func (nthi *NTHi) Hash(canonical bool) <-chan uint64 {
	hashChan := make(chan uint64, bufferSize)
	go func() {
		defer close(hashChan)

		// start the rolling hash
		for {

			// check that rolling can continue
			if nthi.currentIdx >= nthi.maxIdx {
				poolNTHi.Put(nthi)
				return
			}

			// start the hashing
			if nthi.currentIdx != 0 {
				prevBase := (*nthi.seq)[nthi.currentIdx-1]
				endBase := (*nthi.seq)[nthi.currentIdx+nthi.k-1]
				// alg 3. of ntHash paper
				nthi.fh = roL(nthi.fh, 1)
				nthi.fh ^= roL(seedTab[prevBase], nthi.k)
				nthi.fh ^= seedTab[endBase]
				nthi.rh = roR(nthi.rh, 1)
				nthi.rh ^= roR(seedTab[prevBase&offset], 1)
				nthi.rh ^= roL(seedTab[endBase&offset], nthi.k-1)
			}

			// calculate and return the canonical ntHash if requested
			if canonical {
				hashChan <- nthi.getCanonical()
			} else {
				hashChan <- nthi.fh
			}

			// increment the index
			nthi.currentIdx++
		}
	}()
	return hashChan
}

// MultiHash returns a channel to range over the canonical multi ntHash values of a sequence
// canonical is set true to return the canonical k-mers, otherwise the forward hashes are returned
// numMultiHash sets the number of multi hashes to generate for each k-mer
func (nthi *NTHi) MultiHash(canonical bool, numMultiHash uint) <-chan []uint64 {
	hashChan := make(chan []uint64, bufferSize)
	go func() {
		defer close(hashChan)

		// start the rolling hash
		for {

			// check that rolling can continue
			if nthi.currentIdx >= nthi.maxIdx {
				poolNTHi.Put(nthi)
				return
			}

			// start the hashing
			if nthi.currentIdx != 0 {
				prevBase := (*nthi.seq)[nthi.currentIdx-1]
				endBase := (*nthi.seq)[nthi.currentIdx+nthi.k-1]
				// alg 3. of ntHash paper
				nthi.fh = roL(nthi.fh, 1)
				nthi.fh ^= roL(seedTab[prevBase], nthi.k)
				nthi.fh ^= seedTab[endBase]
				nthi.rh = roR(nthi.rh, 1)
				nthi.rh ^= roR(seedTab[prevBase&offset], 1)
				nthi.rh ^= roL(seedTab[endBase&offset], nthi.k-1)
			}

			// set up the return slice
			multiHashes := make([]uint64, numMultiHash)
			if canonical {
				multiHashes[0] = nthi.getCanonical()
			} else {
				multiHashes[0] = nthi.fh
			}

			for i := uint64(1); i < uint64(numMultiHash); i++ {
				tVal := multiHashes[0] * (i ^ uint64(nthi.k)*multiSeed)
				tVal ^= tVal >> multiShift
				multiHashes[i] = tVal
			}

			// send the multihashes for this k-mer
			hashChan <- multiHashes

			// increment the index
			nthi.currentIdx++
		}
	}()
	return hashChan
}

// getCanonical returns the canonical hash value currently held by the iterator
func (nthi *NTHi) getCanonical() uint64 {
	if nthi.rh < nthi.fh {
		return nthi.rh
	}
	return nthi.fh
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
		hv ^= seedTab[seq[i]]
		i++
	}
	return hv
}

// ntr64 generates the ntHash for the reverse strand of the kmer
func ntr64(seq []byte, i, k uint) uint64 {
	var hv uint64
	for i < k {
		hv = roL(hv, 1)
		hv ^= seedTab[seq[k-1-i]&offset]
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
