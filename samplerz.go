package sampler

import (
	"io"
	"math"

	"github.com/holiman/uint256"
)

const (
	// Precision of RCDT
	RCDTprec    uint8 = 72
	RCDTprecLen uint8 = (RCDTprec >> 3)

	inv2sigma2 float64 = 0.15086504887537272 // = 1 / (2 * (math.Pow(1.8205, 2)))

	// ln(2) and 1 / ln(2), with ln the natural logarithm
	LN2  float64 = 0.69314718056
	ILN2 float64 = 1.44269504089
)

// RCDT is the reverse cumulative distribution table of a distribution that
// is very close to a half-Gaussian of parameter MAX_SIGMA.
var RCDT = []*uint256.Int{
	NewBigNumFromHex("0xA3F7F42ED3AC391802"),
	NewBigNumFromHex("0x54D32B181F3F7DDB82"),
	NewBigNumFromHex("0x227DCDD0934829C1FF"),
	NewBigNumFromHex("0xAD1754377C7994AE4"),
	NewBigNumFromHex("0x295846CAEF33F1F6F"),
	NewBigNumFromHex("0x774AC754ED74BD5F"),
	NewBigNumFromHex("0x1024DD542B776AE4"),
	NewBigNumFromHex("0x1A1FFDC65AD63DA"),
	NewBigNumFromHex("0x1F80D88A7B64y28"),
	NewBigNumFromHex("0x1C3FDB2040C69"),
	NewBigNumFromHex("0x12CF24D031FB"),
	NewBigNumFromHex("0x949F8B091F"),
	NewBigNumFromHex("0x3665DA998"),
	NewBigNumFromHex("0xEBF6EBB"),
	NewBigNumFromHex("0x2F5D7E"),
	NewBigNumFromHex("0x7098"),
	NewBigNumFromHex("0xC6"),
	NewBigNumFromHex("0x1"),
}

// C contains the coefficients of a polynomial that approximates exp(-x)
// More precisely, the value:
// (2 ** -63) * sum(C[12 - i] * (x ** i) for i in range(i))
// Should be very close to exp(-x).
// This polynomial is lifted from FACCT: https://doi.org/10.1109/TC.2019.2940949
var C = []*uint256.Int{
	NewBigNumFromInt(0x00000004741183A3),
	NewBigNumFromInt(0x00000036548CFC06),
	NewBigNumFromInt(0x0000024FDCBF140A),
	NewBigNumFromInt(0x0000171D939DE045),
	NewBigNumFromInt(0x0000D00CF58F6F84),
	NewBigNumFromInt(0x000680681CF796E3),
	NewBigNumFromInt(0x002D82D8305B0FEA),
	NewBigNumFromInt(0x011111110E066FD0),
	NewBigNumFromInt(0x0555555555070F00),
	NewBigNumFromInt(0x155555555581FF00),
	NewBigNumFromInt(0x400000000002B400),
	NewBigNumFromInt(0x7FFFFFFFFFFF4800),
	NewBigNumFromInt(0x8000000000000000),
}

type sampler struct {
	y   *uint256.Int
	z   *uint256.Int
	rng io.Reader

	baseSamplerRB []byte // lenght is not checked, but must be RCDTprecLen!
	samplerzRB    []byte // lenght is not checked, but must be 1 byte!
	berexpRB      []byte // lenght is not checked, but must be 1 byte!
}

func newsampler(reader io.Reader) *sampler {
	sp := new(sampler)
	sp.y = new(uint256.Int)
	sp.z = new(uint256.Int)

	sp.rng = reader

	sp.baseSamplerRB = make([]byte, RCDTprecLen)
	sp.samplerzRB = make([]byte, 1)
	sp.berexpRB = make([]byte, 1)

	return sp
}

func (sp *sampler) read(dst []byte) {
	_, err := io.ReadFull(sp.rng, dst)
	if err != nil {
		panic(err)
	}
}

// Require: -
// Ensure: An integer z0 ∈ {0, . . . , 18} such that z ∼ χ ▷ χ is uniquely defined by (3.33)
// 1: u ← UniformBits(72)
// 2: z0 ← 0
// 3: for i = 0, . . . , 17 do
// 4: 	z0 ← z0 + Ju < RCDT[i]K
// 5: return z0
// https://falcon-sign.info/falcon.pdf#57
func (sp *sampler) baseSampler() int {
	var z0 int
	u := sp.y
	sp.read(sp.baseSamplerRB)
	u.SetBytes(sp.baseSamplerRB)
	for _, elt := range RCDT {
		// z0 += 1 if (u < elt)
		if u.Cmp(elt) == -1 {
			z0 += 1
		}
	}
	return z0
}

// Require: Floating-point values x ∈ [0, ln(2)] and ccs ∈ [0, 1]
// Ensure: An integral approximation of 263 · ccs · exp(−x)
// 1: C = [0x00000004741183A3,0x00000036548CFC06,0x0000024FDCBF140A,0x0000171D939DE045,0x0000D00CF58F6F84, 0x000680681CF796E3, 0x002D82D8305B0FEA, 0x011111110E066FD0,0x0555555555070F00, 0x155555555581FF00, 0x400000000002B400, 0x7FFFFFFFFFFF4800,0x8000000000000000]
// 2: y ← C[0]
// 3: z ← ⌊263 · x⌋
// 4: for 1 = 1, . . . , 12 do
// 5: y ← C[u] − (z · y) >> 63
// 6: z ← ⌊263 · ccs⌋
// 7: y ← (z · y) >> 63
// 8: return y
// https://falcon-sign.info/falcon.pdf#d0
func (sp *sampler) approxexp(x, ccs float64) uint64 {
	sp.y.Set(C[0])
	// Since z is positive, int is equivalent to floor
	sp.z.SetUint64(uint64(x * (1 << 63)))
	for _, elt := range C[1:] {
		sp.y.Mul(sp.y, sp.z) // y = z * y
		sp.y.Rsh(sp.y, 63)   // y = y >> 63
		sp.y.Sub(elt, sp.y)  // y = elt - y
	}
	sp.z.SetUint64(uint64(ccs * float64((1<<63)<<1)))
	sp.y.Mul(sp.z, sp.y) // y = z * y
	sp.y.Rsh(sp.y, 63)   // y = y >> 63
	return sp.y.Uint64()
}

// Require: Floating point values x, ccs ≥ 0
// Ensure: A single bit, equal to 1 with probability ≈ ccs · exp(−x)
// 1: s ← ⌊x/ ln(2)⌋
// 2: r ← x − s · ln(2)
// 3: s ← min(s, 63)
// 4: z ← (2 · ApproxExp(r, ccs) − 1) >> s ▷ z ≈ 264−s · ccs · exp(−r) = 264 · ccs · exp(−x)
// 5: i ← 64
// 6: do
// 7: i ← i − 8
// 8: w ← UniformBits(8) − ((z >> i) & 0xFF)
// 9: while ((w = 0) and (i > 0))
// 10: return Jw < 0K ▷ Return 1 with probability 2−64 · z ≈ ccs · exp(−x)
// https://falcon-sign.info/falcon.pdf#cf
func (sp sampler) berexp(x, ccs float64) bool {
	var w int
	s := math.Floor(x * ILN2)
	r := x - s*LN2
	s = Min(s, 63)
	z := (sp.approxexp(r, ccs) - 1) >> int(s)
	for i := 56; i >= -8; i -= 8 {
		sp.read(sp.berexpRB)
		p := int(sp.berexpRB[0])
		w = p - int((z>>uint64(i)))&0xFF
		if w != 0 {
			break
		}
	}
	return w < 0
}

// Given floating-point values mu, sigma (and sigmin),
// output an integer z according to the discrete
// Gaussian distribution D_{Z, mu, sigma}.
//
// Input:
// - the center mu
// - the standard deviation sigma
// - a scaling factor sigmin
// It also takes arguments for underlying functions to prevent unnecessary allocations.
// The inputs MUST verify 1 < sigmin < sigma < MAX_SIGMA.
//
// Output:
// - a sample z from the distribution D_{Z, mu, sigma}.
// https://falcon-sign.info/falcon.pdf#58
func (sp *sampler) Samplerz(mu float64, sigma float64, sigmin float64) int {
	s := int(math.Floor(mu))
	r := mu - float64(s)
	dss := 1 / (2 * sigma * sigma)
	ccs := sigmin / sigma
	for {
		z0 := sp.baseSampler()
		sp.read(sp.samplerzRB)
		b := int(sp.samplerzRB[0])
		b &= 1
		z := float64(b + (2*b-1)*z0)
		x := math.Pow((z-r), 2) * dss
		x -= math.Pow(float64(z0), 2) * inv2sigma2
		if sp.berexp(x, ccs) {
			return s + int(z)
		}
	}
}
