package sampler

import (
	"bytes"
	"encoding/hex"
	"io"

	"github.com/holiman/uint256"
	"golang.org/x/crypto/sha3"
)

func Min(a float64, b float64) float64 {
	if a < b {
		return a
	}
	return b
}

func NewBigNumFromHex(s string) *uint256.Int {
	bn := new(uint256.Int)
	bn.SetFromHex(s)
	return bn
}

func NewBigNumFromInt(i uint64) *uint256.Int {
	bn := uint256.NewInt(i)
	return bn
}

func decodeHexString(hexString string) []byte {
	byteSlice, err := hex.DecodeString(hexString)
	if err != nil {
		panic(err)
	}
	return byteSlice
}

func fromSeedSHAKE(seed []byte) io.Reader {
	shake := sha3.NewShake256()
	_, err := shake.Write(seed)
	if err != nil {
		panic(err) // should never happen
	}
	return shake
}

// Only for testing purposes.
func bytesReader(b []byte) io.Reader {
	return bytes.NewReader(b)
}
