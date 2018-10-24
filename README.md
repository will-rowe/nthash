<div align="center">
  <h1>ntHash</h1>
  <h3><a href="ntHash">ntHash</a> implementation in Go</h3>
  <hr>
  <a href="https://travis-ci.org/will-rowe/ntHash"><img src="https://travis-ci.org/will-rowe/ntHash.svg?branch=master" alt="travis"></a>
  <a href="https://godoc.org/github.com/will-rowe/ntHash"><img src="https://godoc.org/github.com/will-rowe/ntHash?status.svg" alt="GoDoc"></a>
  <a href="https://goreportcard.com/report/github.com/will-rowe/ntHash"><img src="https://goreportcard.com/badge/github.com/will-rowe/ntHash" alt="goreportcard"></a>
  <a href="https://codecov.io/gh/will-rowe/ntHash"><img src="https://codecov.io/gh/will-rowe/ntHash/branch/master/graph/badge.svg" alt="codecov"></a>
</div>

***

## Overview

This is a Go implementation of the [ntHash](https://github.com/bcgsc/ntHash) recursive hash function for hashing all possible k-mers in a DNA/RNA sequence.

For more information, read the ntHash [paper](http://dx.doi.org/10.1093/bioinformatics/btw397) by Mohamadi et al. or check out their C++ [implementation](https://github.com/bcgsc/ntHash).

This implementation was inspired by [Luiz Irber](https://luizirber.org/) and his recent [blog post](https://blog.luizirber.org/2018/09/13/nthash/) on his cool [Rust ntHash implementation](https://github.com/luizirber/nthash).

I have coded this up in Go so that ntHash can be used in my [HULK](https://github.com/will-rowe/hulk) and [GROOT](https://github.com/will-rowe/groot) projects but feel free to use it for yourselves.

## Installation

``` go
go get github.com/will-rowe/ntHash
```

## Example usage

### range over ntHash values for a sequence

``` go
package main

import (
    "log"
    "github.com/will-rowe/ntHash"
)

var (
    sequence = []byte("ACGTCGTCAGTCGATGCAGTACGTCGTCAGTCGATGCAGT")
    kmerSize = 11
)

func main() {
    // create the ntHash iterator using a pointer to the sequence and a k-mer size
    hasher, err := ntHash.New(&sequence, kmerSize)
    // check for errors (e.g. bad k-mer size choice)
    if err != nil {
        log.Fatal(err)
    }
    // collect the hashes by ranging over the hash channel produced by the Hash method
    canonical := true
    for hash := range hasher.Hash(canonical) {
        log.Println(hash)
    }
}
```
