// BSD 3-Clause License
// 
// Copyright (c) 2023, Roy Ward
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// 
// 1. Redistributions of source code must retain the above copyright notice, this
//    list of conditions and the following disclaimer.
// 
// 2. Redistributions in binary form must reproduce the above copyright notice,
//    this list of conditions and the following disclaimer in the documentation
//    and/or other materials provided with the distribution.
// 
// 3. Neither the name of the copyright holder nor the names of its
//    contributors may be used to endorse or promote products derived from
//    this software without specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

#include "pseudo_double.h"
#include <math.h> // Only used to bring in NAN
#include <string.h> // memcpy

pseudo_double_i double_to_pdi(double d) {
	uint64_t i;
	memcpy(&i,&d,sizeof i);   // well-defined in both C and C++, unlike the union
	int32_t raw_exponent=((i>>52)&0x7FF);
	int64_t old_mantissa=(i&0xFFFFFFFFFFFFFULL);
	if(raw_exponent==0x7FF) { // infinity or NaN
#if PD_ERROR_CHECK
		if(old_mantissa==0) {
			PD_DO_ERROR_OVERFLOW;
		} else {
			PD_DO_ERROR_RANGE;
		}
#endif
		return PD_NAN;
	}
	if(raw_exponent==0 && old_mantissa==0) {
		return 0; // +0.0 and -0.0
	}
	// normal: add the implied bit; subnormal: no implied bit, and scale as if raw_exponent were 1
	signed_pd_internal mantissa=(signed_pd_internal)((old_mantissa+(raw_exponent?0x10000000000000ULL:0))<<(PSEUDO_DOUBLE_TOTAL_BITS-54));
	int32_t exponent=(raw_exponent?raw_exponent:1)-0x3FF+2+PSEUDO_DOUBLE_EXP_BIAS;
	if((int64_t)i<0) {
		mantissa=-mantissa;
	}
	// the input is exact, so ties to even; the helper normalises (subnormals, negative powers of two),
	// rounds when the mantissa field is narrower than 53 bits, and handles overflow and underflow
	return pseudo_round_helper(mantissa,exponent,true);
}

// +-m * 2^e, correctly rounded. m must be non-zero.
static inline pseudo_double_i magnitude_scaled_to_pdi(uint64_t m, int negative, int32_t e) {
	int s=PD_CLZ(m);
	m<<=s;                                                     // top bit at 63
	signed_pd_internal vr=(signed_pd_internal)((m>>1)|(m&1));  // positive and normalised; the dropped bit stays as a sticky bit
	if(negative) {
		vr=-vr;
	}
	return pseudo_round_helper(vr,e+65-s+PSEUDO_DOUBLE_EXP_BIAS,true);
}

pseudo_double_i int64_to_pdi(int64_t d) {
	if(d==0) {
		return 0;
	}
	return magnitude_scaled_to_pdi(d<0?0-(uint64_t)d:(uint64_t)d,d<0,0);
}

pseudo_double_i uint64_to_pdi(uint64_t d) {
	if(d==0) {
		return 0;
	}
	return magnitude_scaled_to_pdi(d,0,0);
}

pseudo_double_i int64fixed2_to_pdi(int64_t d, int32_t e) {
	if(d==0) {
		return 0;
	}
	return magnitude_scaled_to_pdi(d<0?0-(uint64_t)d:(uint64_t)d,d<0,e);
}

double pdi_to_double(pseudo_double_i x) {
	if(x==0) {
		return 0.0;
	}
	if(x==PD_NAN) {
		return NAN;
	}
	signed_pd_internal vx=((signed_pd_internal)x)&EXP_MASK_INV;
	uint64_t sgn=(vx<0)?0x8000000000000000ULL:0;
	uint64_t m=(vx<0)?0-(uint64_t)vx:(uint64_t)vx;      // magnitude in [2^62,2^63)
	int s=PD_CLZ(m);
	m<<=s;                                               // top bit at 63: value = (m/2^63)*2^(raw-1023)
	int32_t raw=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS+0x3FF-1-s;
	if(raw>=0x7FF) {                                     // too large for a double: infinity
		uint64_t inf=0x7FF0000000000000ULL|sgn;
		double d;
		memcpy(&d,&inf,sizeof d);
		return d;
	}
	int shift=11;                                        // bits below double's 53-bit significand
	if(raw<1) {                                          // subnormal result: shift further right, encode with raw=1 below
		shift+=1-raw;
		raw=1;
		if(shift>64) {
			return sgn?-0.0:0.0;
		}
	}
	uint64_t sig=(shift==64)?0:(m>>shift);
	uint64_t rem=(shift==64)?m:(m&((1ULL<<shift)-1));
	uint64_t half=1ULL<<(shift-1);
	if(rem>half || (rem==half && (sig&1))) {             // round to nearest, ties to even
		sig++;
	}
	uint64_t bits=((uint64_t)(raw-1)<<52)+sig;           // sig includes the implied bit, so carries move into the exponent
	if(bits>=0x7FF0000000000000ULL) {
		bits=0x7FF0000000000000ULL;                      // overflow: infinity
	}
	bits|=sgn;
	double d;
	memcpy(&d,&bits,sizeof d);
	return d;
}

int64_t pdi_to_int64(pseudo_double_i x) {
	if(x==0) {
		return 0;
	}
	signed_pd_internal vx=((signed_pd_internal)x)&EXP_MASK_INV;
	int32_t shift=PSEUDO_DOUBLE_TOTAL_BITS-((int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS);
	if(shift<0) {                     // |x| >= 2^63: out of range, overflow, or saturate like Rust's `as i64`
#if PD_ERROR_CHECK
		PD_DO_ERROR_INTEGER_OVERFLOW;
#endif
		return (vx<0)?INT64_MIN:INT64_MAX;
	}
	if(shift>=64) {                   // |x| < 1/2
		return 0;
	}
	if(vx>=0) {
		return vx>>shift;
	}
	uint64_t m=0-(uint64_t)vx;        // magnitude, so negatives truncate toward zero like a C cast
	return (int64_t)(0-(m>>shift));
}

int64_t pdi_to_int64fixed2(pseudo_double_i x, int32_t e) {
	if(x==0) {
		return 0;
	}
	if(e<=-1073741824) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_INTEGER_OVERFLOW;
#endif
		return (((signed_pd_internal)x)<0)?INT64_MIN:INT64_MAX;
	}
	if(e>=1073741824) {
		return 0;
	}
	signed_pd_internal vx=((signed_pd_internal)x)&EXP_MASK_INV;
	int32_t shift=PSEUDO_DOUBLE_TOTAL_BITS+e-((int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS);
	if(shift<0) {                     // |x|/2^e >= 2^63: out of range, overflow, or saturate like Rust's `as i64`
#if PD_ERROR_CHECK
		PD_DO_ERROR_INTEGER_OVERFLOW;
#endif
		return (vx<0)?INT64_MIN:INT64_MAX;
	}
	if(shift>=64) {                   // |x|/2^e < 1/2
		return 0;
	}
	if(vx>=0) {
		return vx>>shift;
	}
	uint64_t m=0-(uint64_t)vx;        // magnitude, so negatives truncate toward zero like a C cast
	return (int64_t)(0-(m>>shift));
}

// truncates toward zero; negative values above −1 give 0
uint64_t pdi_to_uint64(pseudo_double_i x) {
	if(x==0) {
		return 0;
	}
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	int32_t exponent=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	if(vx<0) {
#if PD_ERROR_CHECK
		if(exponent>=2 || (exponent==1 && (uint64_t)vx==(1ULL<<63))) { // x <= -1
			PD_DO_ERROR_INTEGER_OVERFLOW;
		}
#endif
		return 0; // -1 < x < 0 truncates to 0; x <= -1 saturates to 0 like Rust's `as u64`
	}
	int32_t shift=PSEUDO_DOUBLE_TOTAL_BITS-exponent;
	if(shift<-1) { // x >= 2^64
#if PD_ERROR_CHECK
		PD_DO_ERROR_INTEGER_OVERFLOW;
#endif
		return UINT64_MAX;
	}
	if(shift==-1) { // 2^63 <= x < 2^64
		return (uint64_t)vx<<1;
	}
	if(shift>=64) { // x < 1/2
		return 0;
	}
	return (uint64_t)vx>>shift;
}

// BEGIN section generated by gen_inv_sqrt_packed.py - do not edit
#define INV_SQRT_C1_SHIFT 44 // c1 is stored as c1*2^20
#define INV_SQRT_C2_SHIFT 41 // c2 is stored as c2*2^23
static const uint64_t inv_sqrt_packed[24]={ // c0:32 | c1:16 | c2:16
	0xFFFEEE8DFECBA550ULL,0xF15B36CCD5BA7D22ULL,0xE4F8AD17B69C6168ULL,0xDA50EC5E9E5D4D94ULL,
	0xD105A5E48B0B3EF8ULL,0xC8D276047B5933F1ULL,0xC18466046E662B71ULL,0xBAF499A9639124C4ULL,
	0xB504D91A5A651F72ULL,0xAF9D3E04528C1B24ULL,0xAAAA993B4BC6179FULL,0xA61D50CE45E114B5ULL,
	0xA1E88EED40B61246ULL,0x9E01A8423C261038ULL,0x9A5FA92938190E7AULL,0x96FAFEAD347B0CFCULL,
	0x93CD33B1313D0BB3ULL,0x90D0BCEE2E510A96ULL,0x8E00D0162BAC099EULL,0x8B594359294508C4ULL,
	0x88D6735C27140804ULL,0x86752E2C2514075BULL,0x8432A216233D06C4ULL,0x820C4F92218D063EULL,
};
// END section generated by gen_inv_sqrt_packed.py - do not edit

// x is a 2.62 unsigned fixed in the range (1,4)
// result is 1.63 unsigned fixed in the range (0.5,1)
uint64_t inv_sqrt64_fixed(uint64_t x) {
	// initial approximation: a quadratic on one of 24 segments of [1,4), good to about 2^-16
	uint64_t w=inv_sqrt_packed[(x>>59)-8];                        // segment from the top 5 bits; c0:32 | c1:16 | c2:16
	uint64_t u=x<<5;                                              // position within the segment, 0.64
	uint64_t c1=((w>>16)&0xFFFF)<<INV_SQRT_C1_SHIFT;              // 0.64
	uint64_t c2=(w&0xFFFF)<<INV_SQRT_C2_SHIFT;                    // 0.64
	uint64_t y=(w&0xFFFFFFFF00000000ULL)-multu64hi(u,c1-multu64hi(u,c2)); // 0.64
	// two Newton-Raphson steps, y=y*(3/2-1/2*x*y*y): error ~2^-16 -> ~2^-31 -> ~2^-60
	y=multu64hi(y,0xC000000000000000ULL-multu64hi(multu64hi(y,y),x))<<1;
	y=multu64hi(y,0xC000000000000000ULL-multu64hi(multu64hi(y,y),x)); // no shift on the last one: result is 1.63
	return y;
}

pseudo_double_i pdi_inv_sqrt(pseudo_double_i x) {
	if(((signed_pd_internal)x)<=0) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_RANGE;
#endif
		return PD_NAN;
	}
	int32_t exponent=x&EXP_MASK;
	uint64_t mantissa=x&EXP_MASK_INV;
	// [01.00 .. 11.11] = [2^0 .. 2^2)  ->  [2^0 .. 2^-1) = [1 .. 0.5)
	if(exponent&1) {
		exponent-=1;
		mantissa<<=1;
	} else {
		if((mantissa<<2)==0) {
			return mantissa+3*(PSEUDO_DOUBLE_EXP_BIAS>>1)+3-(exponent>>1);
		}
	}
	return pseudo_round_helper((signed_pd_internal)inv_sqrt64_fixed(mantissa),3*(PSEUDO_DOUBLE_EXP_BIAS>>1)+2-(exponent>>1),false);
}

pseudo_double_i pdi_sqrt(pseudo_double_i x) {
	if(((signed_pd_internal)x)<0) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_RANGE;
#endif
		return PD_NAN;
	}
	if(x==0) {
		return 0;
	}
	int32_t exponent=x&EXP_MASK;
	uint64_t mantissa=x&EXP_MASK_INV;
	if(exponent&1) {
		exponent-=1;
		mantissa<<=1;
	} else {
		if((mantissa<<2)==0) {
			return mantissa+(PSEUDO_DOUBLE_EXP_BIAS>>1)+1+(exponent>>1);
		}
	}
	// (1,4) * (1,0.5) = (1,2)
	uint64_t y=multu64hi(inv_sqrt64_fixed(mantissa),mantissa)<<1; // 2.62 in (1,2), within a few units of the true root
	int32_t rexp=(PSEUDO_DOUBLE_EXP_BIAS>>1)+1+(exponent>>1);
	// exact correction: q = nearest multiple of the result ulp to y; if y is close to a midpoint,
	// check q against the true root using exact 128-bit squares (a square root can never be an exact tie)
	const uint64_t ulp=1ULL<<PSEUDO_DOUBLE_EXP_BITS, half=ulp>>1;
	uint64_t q=(y+half)&~(ulp-1);
	uint64_t frac=(y+half)&(ulp-1);                              // position of y relative to the rounding boundaries
	if(frac<64 || frac>ulp-64) {                                 // within 64 units of a midpoint: check exactly
		uint64_t t_hi=mantissa>>2, t_lo=mantissa<<62;            // mantissa*2^62: the true root is sqrt(this), in 2.62 units
		uint64_t lo=q-half, hi=q+half;
		uint64_t s_hi=multu64hi(lo,lo), s_lo=lo*lo;
		if(s_hi>t_hi || (s_hi==t_hi && s_lo>t_lo)) {            // lower midpoint is above the root
			q-=ulp;
		} else {
			s_hi=multu64hi(hi,hi);
			s_lo=hi*hi;
			if(s_hi<t_hi || (s_hi==t_hi && s_lo<t_lo)) {        // upper midpoint is below the root
				q+=ulp;
			}
		}
	}
	if(q==(1ULL<<63)) {                                          // rounded up to exactly 2.0
		q=1ULL<<62;
		rexp+=1;
	}
	return q+rexp;
}

// 128-bit signed values are held as (hi,lo): value = hi*2^64 + lo
static inline void pd_add128(int64_t* hi, uint64_t* lo, int64_t bhi, uint64_t blo) {
	uint64_t l=*lo+blo;
	*hi=(int64_t)((uint64_t)*hi+(uint64_t)bhi+(l<blo));
	*lo=l;
}

// BEGIN generated by gen_exp2_table.py - do not edit
static const uint64_t exp2_table[32]={
	9223372036854775808ULL,9425337585438007767ULL,9631725603661576981ULL,9842632931000616231ULL,
	10058158527438640871ULL,10278403519900636096ULL,10503471249702896438ULL,10733467321041881441ULL,
	10968499650544839023ULL,11208678517905446201ULL,11454116617628225966ULL,11704929111906019353ULL,
	11961233684655323370ULL,12223150596734848760ULL,12490802742373206726ULL,12764315706832201104ULL,
	13043817825332782212ULL,13329440243271311081ULL,13621316977754388159ULL,13919584980481119318ULL,
	14224384202002324189ULL,14535857657386837923ULL,14854151493325717732ULL,15179415056705840218ULL,
	15511800964685064948ULL,15851465176301844309ULL,16198567065652879633ULL,16553269496673159336ULL,
	16915738899553466671ULL,17286145348831213003ULL,17664662643191237676ULL,18051468387014017850ULL,
};
// q(r)*2^69 = c1 r + c2 r^2 + ... + c6 r^6
#define EXP2_C1 12786308645202652495ULL
#define EXP2_C2 138481152924917420ULL
#define EXP2_C3 999873132098022ULL
#define EXP2_C4 5414526200327ULL
#define EXP2_C5 23455372821ULL
#define EXP2_C6 85604242ULL
#define EXP2_LOG2E_HI 6653256548922161245LL
#define EXP2_LOG2E_LO 16087983116920987460ULL
#define EXP2_LOG2_10_HI 7659844674706589055LL
#define EXP2_LOG2_10_LO 2636289682573880742ULL
// expm1(x) = x*(1 + u*H(u)), u = 16x; H = h0 + h1 u + ... + h8 u^8 (relative error ~2^-73.2)
#define EXPM1_H0 (4611686018427387904LL)
#define EXPM1_H1 (192153584101141162LL)
#define EXPM1_H2 (6004799503160661LL)
#define EXPM1_H3 (150119987579029LL)
#define EXPM1_H4 (3127499741230LL)
#define EXPM1_H5 (55848209483LL)
#define EXPM1_H6 (872628274LL)
#define EXPM1_H7 (12120806LL)
#define EXPM1_H8 (151508LL)
// END generated by gen_exp2_table.py - do not edit

// 2^x before the final shift: x is 0.64 in [0,1), result is 1.63 in [1,2) (so it needs an unsigned 64-bit value)
static inline uint64_t exp2_64_fixed_inline(uint64_t x) {
	// 2^x = 2^(k/32) * (1+q(r)): table for 2^(k/32), degree 6 polynomial for q, evaluated with Estrin's scheme
	uint64_t t=exp2_table[x>>59];
	uint64_t r=x<<5;
	uint64_t r2=multu64hi(r,r);
	uint64_t a=EXP2_C1+multu64hi(EXP2_C2,r);
	uint64_t b=EXP2_C3+multu64hi(EXP2_C4,r);
	uint64_t c=EXP2_C5+multu64hi(EXP2_C6,r);
	uint64_t r4=multu64hi(r2,r2);
	uint64_t p=a+multu64hi(b,r2)+multu64hi(c,r4);  // q(r)/r * 2^69
	uint64_t q=multu64hi(p,r);                      // q(r) * 2^69
	return t+(multu64hi(t,q)>>5);                   // t + t*q, 1.63
}

// x is a 0.64 unsigned fixed in the range [0,1)
// result is 2.62 unsigned fixed in the range [1,2)
uint64_t exp2_64_fixed(uint64_t x) {
	return exp2_64_fixed_inline(x)>>1;
}


pseudo_double_i pdi_exp2(pseudo_double_i x) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	if(x==0 || e<=-62) {                     // |x| < 2^-63: 2^x rounds to exactly 1
		return (1ULL<<62)+PSEUDO_DOUBLE_EXP_BIAS+2;
	}
	if(e>PSEUDO_DOUBLE_EXP_BITS+1) {         // |x| > 2^EXP_BITS: certain overflow or underflow
		if(vx<0) {
#if PD_ERROR_CHECK
			PD_DO_ERROR_UNDERFLOW;
#endif
			return 0;
		}
#if PD_ERROR_CHECK
		PD_DO_ERROR_OVERFLOW;
#endif
		return PD_NAN;
	}
	// split x = ipart + fraction, with ipart = floor(x) and fraction in [0,1) as 0.64
	int32_t ipart;
	uint64_t fraction;
	if(e>=1) {
		ipart=(int32_t)(vx>>(64-e));
		fraction=(uint64_t)vx<<e;
	} else {
		ipart=(int32_t)(vx>>63);             // 0 or -1
		fraction=(uint64_t)(vx>>-e);
	}
	// exp2_64_fixed gives 2^fraction as 2.62 in [1,2), which is already a normalised mantissa;
	// the helper rounds and handles overflow/underflow of the result exponent
	return pseudo_round_helper((signed_pd_internal)exp2_64_fixed(fraction),ipart+PSEUDO_DOUBLE_EXP_BIAS+2,false);
}
// split t = (h:l)*2^(sh-128), a 128-bit signed value (h:l) with |h:l| in [2^124,2^126], into floor(t) and the fraction
// as 0.64. Returns 0 for that, 1 if |t| < 2^-63 (so 2^t rounds to 1), 2 if t >= 2^EXP_BITS, 3 if t <= -2^EXP_BITS
static inline int exp2_split(int64_t h, uint64_t l, int32_t sh, int32_t* ipart, uint64_t* fraction) {
	if(sh<=-61) {
		return 1;
	}
	if(sh>=PSEUDO_DOUBLE_EXP_BITS+4) {
		return h<0?3:2;
	}
	if(sh>=1) {
		*ipart=(int32_t)(h>>(64-sh));
		*fraction=((uint64_t)h<<sh)|(l>>(64-sh));
	} else {
		*ipart=(int32_t)(h>>63);                                        // 0 or -1
		*fraction=(uint64_t)(h>>-sh);
	}
	return 0;
}

// 2^t, where t = (h:l)*2^(sh-128) as for exp2_split
// shared by pdi_exp, pdi_exp10 and pdi_pow, which all form t in 128 bits so the only rounding is at the end
static inline pseudo_double_i pdi_exp2_wide(int64_t h, uint64_t l, int32_t sh) {
	int32_t ipart;
	uint64_t fraction;
	switch(exp2_split(h,l,sh,&ipart,&fraction)) {
		case 1:
			return (1ULL<<62)+PSEUDO_DOUBLE_EXP_BIAS+2;                 // |t| < 2^-63: 2^t rounds to 1
		case 2:
#if PD_ERROR_CHECK
			PD_DO_ERROR_OVERFLOW;
#endif
			return PD_NAN;
		case 3:
#if PD_ERROR_CHECK
			PD_DO_ERROR_UNDERFLOW;
#endif
			return 0;
	}
	return pseudo_round_helper((signed_pd_internal)exp2_64_fixed(fraction),ipart+PSEUDO_DOUBLE_EXP_BIAS+2,false);
}

// 2^(x*k), with k = (khi:klo)*2^(kexp-128) a 128-bit constant normalised to [2^126,2^127)
static inline pseudo_double_i pdi_exp2_scaled(pseudo_double_i x, int64_t khi, uint64_t klo, int32_t kexp) {
	if(x==0) {
		return (1ULL<<62)+PSEUDO_DOUBLE_EXP_BIAS+2;
	}
	int32_t ex=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;             // x = vx*2^(ex-64)
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	// t = x*k = P*2^(ex+kexp-128), where P (h:l) is the top 128 bits of vx*(khi:klo)
	int64_t h=mults64hi(vx,khi);
	uint64_t l=(uint64_t)vx*(uint64_t)khi;
	int64_t m=mulsu64hi(vx,klo);
	pd_add128(&h,&l,m>>63,(uint64_t)m);
	return pdi_exp2_wide(h,l,ex+kexp);
}

pseudo_double_i pdi_exp(pseudo_double_i x) {
	return pdi_exp2_scaled(x,EXP2_LOG2E_HI,EXP2_LOG2E_LO,2);        // e^x = 2^(x*log2(e))
}

pseudo_double_i pdi_exp10(pseudo_double_i x) {
	return pdi_exp2_scaled(x,EXP2_LOG2_10_HI,EXP2_LOG2_10_LO,3);    // 10^x = 2^(x*log2(10))
}

// e^x - 1. For |x| <= 1/16 a polynomial in relative form keeps full precision as the result gets small; beyond that,
// e^x = 2^t is formed from the exp2 kernel before its final rounding, and 1 is subtracted in fixed point
pseudo_double_i pdi_expm1(pseudo_double_i x) {
	int32_t ex=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;           // x = vx*2^(ex-64)
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	if(ex<=-3) {
		// |x| <= 1/16: expm1(x) = x*(1 + u*H(u)), u = 16x, evaluated with Estrin's scheme
		if(x==0) {
			return 0;
		}
		int64_t u=ex<=-66?(vx>>63):(vx>>(-3-ex));                      // u * 2^63
		int64_t w=mults64hi(u,u);                                      // u^2 * 2^62
		int64_t w2=mults64hi(w,w);                                     // u^4 * 2^60
		int64_t a=EXPM1_H0+mults64hi(EXPM1_H1,u);                      // 2^67
		int64_t b=EXPM1_H2+mults64hi(EXPM1_H3,u);                      // 2^69
		int64_t c=EXPM1_H4+mults64hi(EXPM1_H5,u);                      // 2^71
		int64_t d=EXPM1_H6+mults64hi(EXPM1_H7,u)+mults64hi(EXPM1_H8,w); // 2^73
		int64_t hh=a+mults64hi(b,w)+mults64hi(c+mults64hi(d,w),w2);   // H(u) * 2^67
		int64_t q=(1LL<<62)+(mults64hi(u,hh)>>4);                      // 1 + u*H(u), * 2^62
		return pseudo_round_helper(mults64hi(vx,q),PSEUDO_DOUBLE_EXP_BIAS+ex+2,false);
	}
	// t = x*log2(e) in 128 bits, as in pdi_exp2_scaled
	int64_t h=mults64hi(vx,EXP2_LOG2E_HI);
	uint64_t l=(uint64_t)vx*(uint64_t)EXP2_LOG2E_HI;
	int64_t m=mulsu64hi(vx,EXP2_LOG2E_LO);
	pd_add128(&h,&l,m>>63,(uint64_t)m);
	int32_t ipart;
	uint64_t fraction;
	switch(exp2_split(h,l,ex+2,&ipart,&fraction)) {
		case 2:
#if PD_ERROR_CHECK
			PD_DO_ERROR_OVERFLOW;
#endif
			return PD_NAN;
		case 3:
			return 0xC000000000000000ULL+PSEUDO_DOUBLE_EXP_BIAS+2;      // e^x rounds to 0, so the result is -1
	}                                                                  // (case 1 can't happen for |x| > 1/16)
	uint64_t s=exp2_64_fixed_inline(fraction);                         // 2^fraction, 1.63
	// e^x - 1 = s*2^(ipart-63) - 1; |x| > 1/16 keeps the result away from 0, so absolute precision is enough
	signed_pd_internal vr;
	int32_t exp_max;
	if(ipart>=64) {
		vr=(signed_pd_internal)(s>>1);                                 // the -1 is far below the last bit
		exp_max=ipart+PSEUDO_DOUBLE_EXP_BIAS+2;
	} else if(ipart>=1) {
		vr=(signed_pd_internal)((s-(1ULL<<(63-ipart)))>>1);            // (s*2^ipart - 2^63)*2^-ipart, exact
		exp_max=ipart+PSEUDO_DOUBLE_EXP_BIAS+2;
	} else if(ipart==0) {
		vr=(signed_pd_internal)(s-(1ULL<<63));                         // in (0,1)*2^63
		exp_max=PSEUDO_DOUBLE_EXP_BIAS+1;
	} else if(ipart==-1) {
		vr=(signed_pd_internal)s;                                      // s/2 - 1 = (s - 2^64)*2^-64
		exp_max=PSEUDO_DOUBLE_EXP_BIAS;
	} else {
		vr=(signed_pd_internal)((ipart>-64?s>>-ipart:0)-(1ULL<<63));   // in (-1,-3/4]*2^63
		exp_max=PSEUDO_DOUBLE_EXP_BIAS+1;
	}
	return pseudo_round_helper(vr,exp_max,false);
}

// BEGIN generated by gen_log2_table.py - do not edit
// entry i: { log2(1+i/32) as 0.64, 1/(1+i/32)-1 as signed 0.64 }; entry 32 is m=2 (log2 wraps to 0)
// generated by gen_log2_table.py - do not edit
// entry i: { -log2(c_i) as 0.128 (top 112 bits), with c_i*2^15 (~2^15/(1+i/32)) in the low 16 bits }
static const uint64_t log2_table[33][2]={
	{0ULL,32768ULL},{818952338354772318ULL,8374966259286768671ULL},
	{1613810734155728114ULL,10207760053733259384ULL},{2385129909880439749ULL,1275646822876476679ULL},
	{3134664534162957034ULL,4553028396277199303ULL},{3863631110689725337ULL,17248725768077012660ULL},
	{4573556910156900338ULL,5822213613265710026ULL},{5264310230205635445ULL,16783178984127817991ULL},
	{5938931262175287918ULL,16574822553249998438ULL},{6595695702403272484ULL,1044430739071460327ULL},
	{7237080614488995453ULL,14550440788060823942ULL},{7863729826320181908ULL,1386051784181636929ULL},
	{8475322052029856596ULL,2718517490555968791ULL},{9072732869993225441ULL,266594463135062790ULL},
	{9658165691507551909ULL,15430077011846584587ULL},{10230510894422911743ULL,7783732352430135078ULL},
	{10791059629171537429ULL,7016836197740008789ULL},{11338786057124349886ULL,6981792152318333848ULL},
	{11876441236189325642ULL,6102940079770980844ULL},{12404464277676035218ULL,3293347085104730192ULL},
	{12920730699165356579ULL,11816404982556675781ULL},{13428371462146184988ULL,16553329721659837768ULL},
	{13925318077683264137ULL,16085359501358222298ULL},{14413568045068074423ULL,17932952614796741241ULL},
	{14892460506987555979ULL,3898792729163155749ULL},{15364210453677207441ULL,17354544718346733532ULL},
	{15826803556703254300ULL,4444950315043538591ULL},{16282601093972579564ULL,456271980907677036ULL},
	{16729584805695595022ULL,9660409584621470788ULL},{17168718146006484123ULL,15619320035400434470ULL},
	{17601053739889444981ULL,959456717530743313ULL},{18027734158009302557ULL,7636027819432296708ULL},
	{0ULL,16384ULL},
};
// log2(1+r)/r = b0 + u*(b1 + b2 u + ... + b8 u^7), u = r*2^6; interpolation error ~2^-65.3
#define LOG2_C0 (6653256548922161246LL)
#define LOG2_C1 (-3326628274461080540LL)
#define LOG2_C2 (69304755717939177LL)
#define LOG2_C3 (-1624330212143633LL)
#define LOG2_C4 (40608255303606LL)
#define LOG2_C5 (-1057506584682LL)
#define LOG2_C6 (28326069129LL)
#define LOG2_C7 (-774881481LL)
#define LOG2_C8 (21524701LL)
#define LOG2_C0LO (-2358760956788564156LL) // b0*2^126 - LOG2_C0*2^64
// ln(2) and 2*log10(2) as 1.63
#define LOG2_TO_LN (6393154322601327830LL)
#define LOG2_TO_LOG10_X2 (5553023288523357132LL)
// END generated by gen_log2_table.py - do not edit

// r is u=r*2^6 as signed 1.63 (|r| < 2^-6), r2 = mults64hi(r,r) = u^2 * 2^62
// result is the tail of log2(1+r)/r = b0 + u*(b1 + b2 u + ... + b8 u^7), i.e. (b1 + ... + b8 u^7) * 2^(68-s)
// s is a constant: 5 for log2_64_fixed (enough precision, and no shift needed after), 0 for pdi_pow
static inline int64_t log2_tail(int64_t r, int64_t r2, int s) {
	int64_t r4=mults64hi(r2,r2);                                  // u^4 * 2^60
	int64_t a=(LOG2_C1>>s)+mults64hi(LOG2_C2>>s,r);               // 2^(68-s)
	int64_t b=(LOG2_C3>>s)+mults64hi(LOG2_C4>>s,r);               // 2^(70-s)
	int64_t c=(LOG2_C5>>s)+mults64hi(LOG2_C6>>s,r);               // 2^(72-s)
	int64_t d=(LOG2_C7>>s)+mults64hi(LOG2_C8>>s,r);               // 2^(74-s)
	return a+mults64hi(b,r2)+mults64hi(c+mults64hi(d,r2),r4);     // Estrin's scheme
}

// r = (1+x)*c - 1 scaled by 2^69, where c = cs/2^15 is the short reciprocal from the table
// exact when the low 10 bits of x are 0 (always true for a pdi mantissa); |r| < 2^-6 so the bits above 2^63 drop out
static inline int64_t log2_reduce(uint64_t x, uint64_t cs) {
	return (int64_t)((x>>10)*cs+(((x&1023)*cs)>>10)+((cs-32768)<<54)); // (c-1)*2^69 + x*c*2^5
}

static inline uint64_t log2_64_fixed_inline(uint64_t x) {
	// log2(1+x) = -log2(c_i) + log2(1+r): c_i ~ 1/(1+i/32), with 1+i/32 the nearest 32nd to 1+x, and r = (1+x)*c_i - 1
	const uint64_t* w=log2_table[((x>>58)+1)>>1];                 // i = 32 when 1+x rounds up to 2
	int64_t r=log2_reduce(x,w[1]&0xFFFF);
	int64_t r2=mults64hi(r,r);
	int64_t y=(mults64hi(r,LOG2_C0)+4)+mults64hi(r2,log2_tail(r,r2,5)); // r*P(r) = r*b0 + r*u*tail, * 2^67 (+4 to round)
	return w[0]+(uint64_t)(y>>3);
}

// x is a 0.64 unsigned fixed in the range [0,1)
// result is log2(1+x) as 0.64 unsigned fixed in the range [0,1)
uint64_t log2_64_fixed(uint64_t x) {
	return log2_64_fixed_inline(x);
}

// log2(1+r) for |r| <= 2^-6, before rounding (vr*2^(exp_max-PSEUDO_DOUBLE_EXP_BIAS-64), the form pseudo_round_helper
// takes). r = rn*2^-(69+s) with rn normalised, so the result carries full relative precision; rr = r*2^69 is only used
// in the tail, where absolute precision is enough. log2(1+r) = r*b0 + r*u*tail
static inline signed_pd_internal log2_small_u(int64_t rn, int32_t s, int64_t rr, int32_t* exp_max) {
	*exp_max=PSEUDO_DOUBLE_EXP_BIAS-3-s;
	return mults64hi(rn,LOG2_C0)+mults64hi(mults64hi(rn,rr),log2_tail(rr,mults64hi(rr,rr),5));
}

// e + log2(1+f) before rounding, f as 0.64: puts the integer e above the 0.64 fraction, using as few bits for e as
// possible. Absolute precision, so not for results near 0 (use log2_small_u)
static inline signed_pd_internal log2_general_u(int32_t e, uint64_t f, int32_t* exp_max) {
	uint64_t lo=log2_64_fixed_inline(f);
	uint32_t a=(uint32_t)(e<0?~e:e);
	int32_t nb=63-PD_CLZ(((uint64_t)a<<1)|1);                       // bits needed for e, excluding the sign bit
	if(nb==0 && (int64_t)(lo^(uint64_t)(int64_t)e)>=0) {
		*exp_max=PSEUDO_DOUBLE_EXP_BIAS;
		return (signed_pd_internal)lo;                              // e is 0 or -1 and |result| < 1/2: lo as signed 0.64 is exactly the result
	}
	*exp_max=PSEUDO_DOUBLE_EXP_BIAS+1+nb;
	return (signed_pd_internal)(((uint64_t)(int64_t)e<<(63-nb))|(lo>>(nb+1)));
}

// log2(x) for x > 0, before rounding: the result is vr*2^(exp_max-PSEUDO_DOUBLE_EXP_BIAS-64), the form pseudo_round_helper takes
// vr is not necessarily normalised, and is 0 for x = 1
static inline int32_t pdi_log2_unrounded(pseudo_double_i x, signed_pd_internal* vr_out) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS-2;      // x = (1+f)*2^e
	uint64_t f=(uint64_t)(x&EXP_MASK_INV)<<2;                       // f as 0.64
	int32_t exp_max=0;
	if((e==0 && f<(1ULL<<58)) || (e==-1 && f>=0-(1ULL<<58))) {
		// x in [1-2^-7,1+2^-6): the result is small, so it needs relative precision
		int64_t r=(int64_t)(e==0?f<<5:f<<4);                        // r = x-1, scaled by 2^69 (exact)
		if(r==0) {
			*vr_out=0;                                              // log(1) = 0
			return 0;
		}
		int32_t s=PD_CLZ(r>0?r:~r)-1;
		*vr_out=log2_small_u((int64_t)((uint64_t)r<<s),s,r,&exp_max);
	} else {
		*vr_out=log2_general_u(e,f,&exp_max);
	}
	return exp_max;
}

// log2(x) * k * 2^(kexp-63), with a single rounding at the end; k=0 means log2(x) itself
static inline pseudo_double_i pdi_log2_scaled(pseudo_double_i x, int64_t k, int32_t kexp) {
	if((signed_pd_internal)x<=0) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_RANGE;
#endif
		return PD_NAN;
	}
	signed_pd_internal vr;
	int32_t exp_max=pdi_log2_unrounded(x,&vr);
	if(vr==0) {
		return 0;                                                   // log(1) = 0
	}
	if(k!=0) {
		int32_t z=PD_CLZ(vr>0?vr:~vr)-1;                           // normalise first, so the product keeps full precision
		vr=mults64hi((signed_pd_internal)((uint64_t)vr<<z),k);
		exp_max+=1+kexp-z;
	}
	return pseudo_round_helper(vr,exp_max,false);
}

pseudo_double_i pdi_log2(pseudo_double_i x) {
	return pdi_log2_scaled(x,0,0);
}

pseudo_double_i pdi_log(pseudo_double_i x) {
	return pdi_log2_scaled(x,LOG2_TO_LN,0);          // ln(x) = log2(x)*ln(2)
}

pseudo_double_i pdi_log10(pseudo_double_i x) {
	return pdi_log2_scaled(x,LOG2_TO_LOG10_X2,-1);   // log10(x) = log2(x)*log10(2)
}

// extended precision log2, for pdi_pow
// log2(x) for x > 0, to about 2^-70 relative: the result is (hi:lo)*2^(exp-128), with |hi:lo| in [2^126,2^127]
// returns exp, or INT32_MIN for x = 1 (log2 = 0)
static inline int32_t pdi_log2_extended(pseudo_double_i x, int64_t* hi_out, uint64_t* lo_out) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS-2;      // x = (1+f)*2^e
	uint64_t f=(uint64_t)(x&EXP_MASK_INV)<<2;                       // f as 0.64
	int64_t hi;
	uint64_t lo;
	int32_t ex;
	if((e==0 && f<(1ULL<<58)) || (e==-1 && f>=0-(1ULL<<58))) {
		// x in [1-2^-7,1+2^-6): log2(x) = r*b0 + r*u*tail with r = x-1 normalised, all scaled by 2^(131+s)
		int64_t r=(int64_t)(e==0?f<<5:f<<4);                        // r = x-1, scaled by 2^69 (exact)
		if(r==0) {
			return INT32_MIN;
		}
		int32_t s=PD_CLZ(r>0?r:~r)-1;
		int64_t rn=(int64_t)((uint64_t)r<<s);
		hi=mults64hi(rn,LOG2_C0);                                   // r*b0, with b0 to 128 bits
		lo=(uint64_t)rn*(uint64_t)LOG2_C0;
		int64_t m=mults64hi(rn,LOG2_C0LO);
		pd_add128(&hi,&lo,m>>63,(uint64_t)m);
		int64_t t=mults64hi(mults64hi(rn,r),log2_tail(r,mults64hi(r,r),0)); // r*u*tail * 2^(72+s)
		pd_add128(&hi,&lo,t>>5,(uint64_t)t<<59);
		ex=-3-s;
	} else {
		// log2(1+f) = -log2(c_i) + log2(1+r) as 0.128, then e is put above it
		const uint64_t* w=log2_table[((f>>58)+1)>>1];
		int64_t r=log2_reduce(f,w[1]&0xFFFF);                       // exact
		int64_t r2=mults64hi(r,r);
		int64_t yh=mults64hi(r,LOG2_C0);                            // r*b0 * 2^131, with b0 to 128 bits
		uint64_t yl=(uint64_t)r*(uint64_t)LOG2_C0;
		int64_t m=mults64hi(r,LOG2_C0LO);
		pd_add128(&yh,&yl,m>>63,(uint64_t)m);
		yl=(yl>>3)|((uint64_t)yh<<61);                              // 2^131 -> 2^128
		yh>>=3;
		int64_t t=mults64hi(r2,log2_tail(r,r2,0));                  // r*u*tail * 2^72
		pd_add128(&yh,&yl,t>>8,(uint64_t)t<<56);
		uint64_t fl=(w[1]&~0xFFFFULL)+yl;                           // -log2(c_i) + log2(1+r), mod 1
		uint64_t fh=w[0]+(uint64_t)yh+(fl<yl);
		uint32_t a=(uint32_t)(e<0?~e:e);
		int32_t nb=63-PD_CLZ(((uint64_t)a<<1)|1);                   // bits needed for e, excluding the sign bit
		if(nb==0 && (int64_t)(fh^(uint64_t)(int64_t)e)>=0) {
			hi=(int64_t)fh;                                         // e is 0 or -1 and |result| < 1/2
			lo=fl;
			ex=0;
		} else {
			hi=(int64_t)(((uint64_t)(int64_t)e<<(63-nb))|(fh>>(nb+1)));
			lo=(fh<<(63-nb))|(fl>>(nb+1));
			ex=nb+1;
		}
	}
	int32_t z=PD_CLZ(hi>0?hi:~hi)-1;                                // normalise
	if(z>0) {
		hi=(int64_t)(((uint64_t)hi<<z)|(lo>>(64-z)));
		lo<<=z;
		ex-=z;
	}
	*hi_out=hi;
	*lo_out=lo;
	return ex;
}

// ln(1+x), accurate for small x as well as large
pseudo_double_i pdi_ln_1p(pseudo_double_i x) {
	int32_t ex=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;           // x = vx*2^(ex-64)
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	signed_pd_internal vr;
	int32_t exp_max;
	if(ex<=-5) {
		// |x| <= 2^-6: log2(1+x) directly from r = x, so 1+x is never formed
		if(x==0) {
			return 0;
		}
		vr=log2_small_u(vx,-5-ex,ex<=-69?(vx>>63):(vx>>(-5-ex)),&exp_max);
	} else {
		if(vx<0 && ex>=2) {
#if PD_ERROR_CHECK
			PD_DO_ERROR_RANGE;                                         // x <= -1
#endif
			return PD_NAN;
		}
		int32_t e;
		uint64_t f;
		if(ex>=64) {
			e=ex-2;                                                    // x >= 2^62: the 1 is far below the last bit
			f=(uint64_t)vx<<2;
		} else {
			// 1+x exactly, as hi:lo with 64 fraction bits: 2^64 + vx*2^ex (exact: the low bits of vx are 0)
			int64_t hi;
			uint64_t lo;
			if(ex>=1) {
				hi=vx>>(64-ex);
				lo=(uint64_t)vx<<ex;
			} else {
				hi=vx>>63;
				lo=(uint64_t)(vx>>-ex);
			}
			hi+=1;                                                     // positive, since x > -1
			// split into 2^e*(1+f)
			if(hi==0) {
				int32_t c=PD_CLZ(lo);
				e=-1-c;
				f=c==63?0:lo<<(c+1);
			} else {
				int32_t c=PD_CLZ((uint64_t)hi);
				e=63-c;
				f=c==63?lo:((uint64_t)hi<<(c+1))|(lo>>(63-c));
			}
		}
		vr=log2_general_u(e,f,&exp_max);
	}
	int32_t z=PD_CLZ(vr>0?vr:~vr)-1;                                   // normalise, then times ln(2), then round once
	vr=mults64hi((signed_pd_internal)((uint64_t)vr<<z),LOG2_TO_LN);
	return pseudo_round_helper(vr,exp_max+1-z,false);
}

// x^y = 2^(y*log2(x)): log2(x) is computed to about 2^-70 relative (pdi_log2_extended) and
// multiplied by y in 128 bits, because any error in t = y*log2(x) is magnified by |y|
pseudo_double_i pdi_pow(pseudo_double_i x, pseudo_double_i y) {
	const pseudo_double_i one=(1ULL<<62)+PSEUDO_DOUBLE_EXP_BIAS+2;
	if(y==0) {
		return one;                                                     // x^0 = 1, including 0^0
	}
	int32_t ey=(int32_t)(y&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;             // y = vy*2^(ey-64)
	signed_pd_internal vy=(signed_pd_internal)(y&EXP_MASK_INV);
	bool negate=false;
	if((signed_pd_internal)x<=0) {
		if(x==0) {
			if(vy>0) {
				return 0;                                               // 0^y = 0 for y > 0
			}
#if PD_ERROR_CHECK
			PD_DO_ERROR_RANGE;                                          // 0^y is infinite for y < 0
#endif
			return PD_NAN;
		}
		// negative x is only defined for integer y: |x|^y, negated if y is odd
		// for 1 <= ey <= 63 the bits of vy below the units bit are vy<<ey; for ey >= 64, y is an even integer
		if(ey<1 || (ey<64 && ((uint64_t)vy<<ey)!=0)) {
#if PD_ERROR_CHECK
			PD_DO_ERROR_RANGE;
#endif
			return PD_NAN;
		}
		negate=(ey<64) && ((vy>>(64-ey))&1);
		x=pdi_neg(x);
	}
	int64_t lh;
	uint64_t ll;
	int32_t el=pdi_log2_extended(x,&lh,&ll);                            // log2(x) = (lh:ll)*2^(el-128)
	if(el==INT32_MIN) {
		return negate?pdi_neg(one):one;                                 // x = 1
	}
	// t = y*log2(x) = P*2^(sh-128), where P (h:l) is the top 128 bits of vy*(lh:ll); |P| is in [2^124,2^126]
	int32_t sh=el+ey;
	signed_pd_internal h=mults64hi(vy,lh);
	uint64_t l=(uint64_t)vy*(uint64_t)lh;
	int64_t m=mulsu64hi(vy,ll);
	pd_add128(&h,&l,m>>63,(uint64_t)m);
	pseudo_double_i r=pdi_exp2_wide(h,l,sh);
	return negate?pdi_neg(r):r;
}

// #if __GNUG__
// // static const pseudo_double_i pdi_inv_log_2_e =int64fixed10_to_pdi(6931471805599453094,-19);
// // static const pseudo_double_i pdi_inv_log_2_10=int64fixed10_to_pdi(3010299956639811952,-19);
// static const pseudo_double_i pdi_1_div_tau=pdi_div(uint64_to_pdi(1),int64fixed10_to_pdi(6283185307179586477,-18));
// static const pseudo_double_i pdi_tau=int64fixed10_to_pdi(6283185307179586477,-18);
// #endif
//

// BEGIN generated by gen_sin_table.py - do not edit
// entry k: { sin(k/32 quarter turn), cos(k/32 quarter turn) } as 1.63
static const uint64_t sin_table[33][2]={
	{0ULL,9223372036854775808ULL},{452569415305004282ULL,9212262081300395919ULL},
	{904048551248139761ULL,9178958979493303928ULL},{1353349755046017570ULL,9123542961522760326ULL},
	{1799390620744551094ULL,9046147529429926059ULL},{2241096596829706927ULL,8946959135589751977ULL},
	{2677403574916221777ULL,8826216733530876278ULL},{3107260453277907452ULL,8684211202275644159ULL},
	{3529631669043774883ULL,8521284645587064995ULL},{3943499692959694935ULL,8337829567810882499ULL},
	{4347867480705496636ULL,8134287928298226505ULL},{4741760874862065242ULL,7911150076686825028ULL},
	{5124230951741890994ULL,7668953571605777421ULL},{5494356307429348228ULL,7408281885649732303ULL},
	{5851245277523433568ULL,7129762999742300884ULL},{6194038085235408523ULL,6834067890275007352ULL},
	{6521908912666391106ULL,6521908912666391106ULL},{6834067890275007352ULL,6194038085235408523ULL},
	{7129762999742300884ULL,5851245277523433568ULL},{7408281885649732303ULL,5494356307429348228ULL},
	{7668953571605777421ULL,5124230951741890994ULL},{7911150076686825028ULL,4741760874862065242ULL},
	{8134287928298226505ULL,4347867480705496636ULL},{8337829567810882499ULL,3943499692959694935ULL},
	{8521284645587064995ULL,3529631669043774883ULL},{8684211202275644159ULL,3107260453277907452ULL},
	{8826216733530876278ULL,2677403574916221777ULL},{8946959135589751977ULL,2241096596829706927ULL},
	{9046147529429926059ULL,1799390620744551094ULL},{9123542961522760326ULL,1353349755046017570ULL},
	{9178958979493303928ULL,904048551248139761ULL},{9212262081300395919ULL,452569415305004282ULL},
	{9223372036854775808ULL,0ULL},
};
// sin(t)/u = a0 + a1 z + a2 z^2 + a3 z^3 (abs error ~2^-73.6), (1-cos(t))/u^2 = e0 + ... (abs error ~2^-82.3)
#define SIN_A0 (1811004864519280711LL)
#define SIN_A1 (-727290917205962LL)
#define SIN_A2 (87622969193LL)
#define SIN_A3 (-5026911LL)
#define SIN_E0 (44448746703040243LL)
#define SIN_E1 (-8925202353581LL)
#define SIN_E2 (716863740LL)
#define SIN_E3 (-30845LL)
#define SIN_INV_TAU_0 (5871781006564002452ULL)
#define SIN_INV_TAU_1 (18308165927316385504ULL)
#define SIN_INV_TAU_2 (7904181063698728993ULL)
#define SIN_2_DIV_PI (5871781006564002453LL) // 2/pi as 1.63, for small x: u = x*128/pi
// END generated by gen_sin_table.py - do not edit

// b is u as signed 1.63 (u in [-1,1] is the offset in 1/64ths of a quarter turn, so the angle is u*pi/128),
// z = mults64hi(b,b) = u^2 * 2^62
// result is sin(angle)/u scaled by 2^66
static inline int64_t sin_poly(int64_t z, int64_t z2) {
	return (SIN_A0+mults64hi(SIN_A1,z))+mults64hi(SIN_A2+mults64hi(SIN_A3,z),z2);
}

// x:xlo is 2.126 in [0,1] quarter turns (xlo is 0 except for the radian functions, where the turn fraction is not exact).
// Returns sin(a_k) and cos(a_k) (1.63) for a_k the nearest 32nd of a quarter turn,
// and corrections so that sin(x quarter turns) = (2*s + ds) * 2^-64 and cos(x quarter turns) = (2*c + dc) * 2^-64
// (the extra bit matters to the pdi functions). Callers that only need one of sin/cos let the compiler drop the other.
static inline void sin_cos_rev_kernel(uint64_t x, uint64_t xlo, uint64_t* s, int64_t* ds, uint64_t* c, int64_t* dc) {
	// sin(a_k + t) = sin(a_k) + cos(a_k)*sin(t) - sin(a_k)*(1-cos(t))
	// cos(a_k + t) = cos(a_k) - sin(a_k)*sin(t) - cos(a_k)*(1-cos(t))
	const uint64_t* w=sin_table[((x>>56)+1)>>1];
	int64_t b=(int64_t)((x<<7)|(xlo>>57));                             // u*2^63: |t| <= 1/64 so the bits above drop out
	int64_t z=mults64hi(b,b);                                          // u^2 * 2^62
	int64_t z2=mults64hi(z,z);                                         // u^4 * 2^60
	int64_t sp=sin_poly(z,z2);                                         // sin(t)/u * 2^66
	int64_t ep=(SIN_E0+mults64hi(SIN_E1,z))+mults64hi(SIN_E2+mults64hi(SIN_E3,z),z2); // (1-cos(t))/u^2 * 2^67
	int64_t cu=mulsu64hi(b,w[1]);                                     // cos(a_k)*u * 2^62
	int64_t su=mulsu64hi(b,w[0]);                                     // sin(a_k)*u * 2^62
	int64_t sz=(int64_t)multu64hi(w[0],(uint64_t)z);                   // sin(a_k)*u^2 * 2^61
	int64_t cz=(int64_t)multu64hi(w[1],(uint64_t)z);                   // cos(a_k)*u^2 * 2^61
	*s=w[0];
	*ds=mults64hi(cu,sp)-mults64hi(sz,ep);                             // * 2^64
	*c=w[1];
	*dc=-mults64hi(su,sp)-mults64hi(cz,ep)-1;                          // * 2^64 (-1 recentres the two negated truncations)
}

// x is a 2.62 unsigned fixed in the range [0,1]
// calculate sin_rev_64_fixed(x) = sin(x quarter turns)
// result is 2.62 unsigned fixed in the range [0,1]
uint64_t sin_rev_64_fixed(uint64_t x) {
	uint64_t s,c;
	int64_t ds,dc;
	sin_cos_rev_kernel(x,0,&s,&ds,&c,&dc);
	return (s+(uint64_t)((ds+2)>>1)+1)>>1;                             // s + ds/2 rounded, then 1.63 -> 2.62 rounded
}

// as sin_rev_64_fixed, but both sin and cos of x quarter turns
void sin_cos_rev_64_fixed(uint64_t x, uint64_t* sin_out, uint64_t* cos_out) {
	uint64_t s,c;
	int64_t ds,dc;
	sin_cos_rev_kernel(x,0,&s,&ds,&c,&dc);
	*sin_out=(s+(uint64_t)((ds+2)>>1)+1)>>1;
	*cos_out=(c+(uint64_t)((dc+2)>>1)+1)>>1;
}

// |sin(x turns)| for |x| <= 2^-8 turn, before rounding (see pseudo_round_helper); ax = |vx|.
// Worked from the mantissa directly, since x as a fixed point turn would lose relative precision.
// With u = x*2^8 (so x turns = u*pi/128 radians), sin = u*sin_poly: ax*sin_poly*2^(e-58)
static inline signed_pd_internal sin_rev_tiny_u(int32_t e, uint64_t ax, int32_t* exp_max) {
	uint64_t b=e<=-70?0:(ax>>(-7-e));                                  // |u| * 2^63
	int64_t z=(int64_t)multu64hi(b,b);
	*exp_max=PSEUDO_DOUBLE_EXP_BIAS+e+6;
	return (signed_pd_internal)multu64hi(ax,(uint64_t)sin_poly(z,mults64hi(z,z)));
}

// sin(x turns) for |x| <= 2^-8 turn: the magnitude is rounded, so sin(-x) = -sin(x) exactly
static inline pseudo_double_i sin_rev_tiny(pseudo_double_i x, int32_t e, signed_pd_internal vx) {
	int32_t exp_max;
	signed_pd_internal vr=sin_rev_tiny_u(e,vx<0?0-(uint64_t)vx:(uint64_t)vx,&exp_max);
	pseudo_double_i r=pseudo_round_helper(vr,exp_max,false);
	return vx<0?pdi_neg(r):r;
}

// x mod 1 turn, folded into the first quarter turn: returns t in [0,1] quarter turns as 2.62, with
// sin(x turns) = (*negative ? -1 : 1) * sin(t quarter turns) and cos(x turns) = (*reflected ? -1 : 1) * cos(t quarter turns)
static inline uint64_t sin_rev_fold(int32_t e, signed_pd_internal vx, bool* negative, bool* reflected) {
	uint64_t t;                                                        // x mod 1 turn, as 0.64 (exact unless x is tiny)
	if(e>=64) {
		t=0;                                                           // a whole number of turns
	} else if(e>=1) {
		t=(uint64_t)vx<<e;
	} else if(e>=-63) {
		t=(uint64_t)(vx>>-e);
	} else {
		t=(uint64_t)(vx>>63);
	}
	*negative=(int64_t)t<0;                                            // second half turn: sin(x) = -sin(-x), cos(x) = cos(-x)
	if(*negative) {
		t=0-t;
	}
	*reflected=t>(1ULL<<62);                                           // second quarter: sin(x) = sin(1/2-x), cos(x) = -cos(1/2-x)
	if(*reflected) {
		t=(1ULL<<63)-t;
	}
	return t;
}

// sin(t:tlo quarter turns) for t:tlo (2.126) below 2^56 (within 1/64 quarter turn of 0), before rounding; normalising u
// keeps full relative precision. tlo is 0 except for the radian functions. Returns 0 for t:tlo = 0
static inline signed_pd_internal sin_rev_small128_u(uint64_t t, uint64_t tlo, int32_t* exp_max) {
	uint64_t h=(t<<7)|(tlo>>57);                                       // u * 2^63 (u in [0,1)), with the next 64 bits in l
	uint64_t l=tlo<<7;
	*exp_max=0;
	if(h==0) {
		return 0;                                                      // only for t = 0 (u below 2^-63 can't occur)
	}
	int32_t s=PD_CLZ(h)-1;
	uint64_t hn=s==0?h:(h<<s)|(l>>(64-s));
	int64_t z=mults64hi((int64_t)h,(int64_t)h);
	*exp_max=PSEUDO_DOUBLE_EXP_BIAS-1-s;
	return mults64hi((int64_t)hn,sin_poly(z,mults64hi(z,z)));         // * 2^(65+s)
}

static inline pseudo_double_i sin_rev_small128(uint64_t t, uint64_t tlo) {
	int32_t exp_max;
	signed_pd_internal vr=sin_rev_small128_u(t,tlo,&exp_max);
	return pseudo_round_helper(vr,exp_max,false);
}

static inline pseudo_double_i sin_rev_small(uint64_t t) {
	return sin_rev_small128(t,0);
}

// r is a 0.64 value in [sin(pi/128),1] where exactly 1 has wrapped to 0: returns it before rounding
static inline signed_pd_internal sin_rev_round_u(uint64_t r, int32_t* exp_max) {
	if(r==0) {
		*exp_max=PSEUDO_DOUBLE_EXP_BIAS+2;
		return 1LL<<62;                                                // 1
	}
	if((int64_t)r>0) {
		*exp_max=PSEUDO_DOUBLE_EXP_BIAS;
		return (signed_pd_internal)r;                                  // < 1/2: all 64 bits fit
	}
	*exp_max=PSEUDO_DOUBLE_EXP_BIAS+1;
	return (signed_pd_internal)(r>>1);
}

static inline pseudo_double_i sin_rev_round(uint64_t r) {
	int32_t exp_max;
	signed_pd_internal vr=sin_rev_round_u(r,&exp_max);
	return pseudo_round_helper(vr,exp_max,false);
}

// tan = (vs*2^es)/(vc*2^ec), from the unrounded magnitudes of sin and cos, with a single rounding
static inline pseudo_double_i tan_rev_parts(signed_pd_internal vs, int32_t es, signed_pd_internal vc, int32_t ec, bool negative) {
	if(vc==0) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_RANGE;                                             // an exact pole: an odd number of quarter turns
#endif
		return PD_NAN;
	}
	if(vs==0) {
		return 0;
	}
	int32_t zs=PD_CLZ(vs)-1;
	int32_t zc=PD_CLZ(vc)-1;
	pseudo_double_i r=pseudo_round_helper(divs64hi((vs<<zs)>>2,vc<<zc),es-zs-ec+zc+PSEUDO_DOUBLE_EXP_BIAS+2,false);
	return negative?pdi_neg(r):r;
}

// each result is rounded as a magnitude and then given its sign, so sin(-x) = -sin(x), cos(-x) = cos(x),
// sin(1/2-x) = sin(x) and cos(1/2-x) = -cos(x) exactly; pdi_sin_cos_rev gives the same bits as pdi_sin_rev and pdi_cos_rev
#define PD_SIN_SMALL (1ULL<<56)                                        // below this (in quarter turns) a result needs relative precision
#define PD_QUARTER (1ULL<<62)

pseudo_double_i pdi_sin_rev(pseudo_double_i x) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;           // x = vx*2^(e-64)
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	if(e<=-7) {
		return sin_rev_tiny(x,e,vx);
	}
	bool negative,reflected;
	uint64_t t=sin_rev_fold(e,vx,&negative,&reflected);
	pseudo_double_i r;
	if(t<PD_SIN_SMALL) {
		r=sin_rev_small(t);
	} else {
		uint64_t s,c;
		int64_t ds,dc;
		sin_cos_rev_kernel(t,0,&s,&ds,&c,&dc);
		r=sin_rev_round(2*s+(uint64_t)ds);
	}
	return negative?pdi_neg(r):r;
}

pseudo_double_i pdi_cos_rev(pseudo_double_i x) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	bool negative,reflected;
	uint64_t t=sin_rev_fold(e,vx,&negative,&reflected);
	pseudo_double_i r;
	if(t>PD_QUARTER-PD_SIN_SMALL) {
		r=sin_rev_small(PD_QUARTER-t);                                 // cos(t) = sin(1/4-t)
	} else {
		uint64_t s,c;
		int64_t ds,dc;
		sin_cos_rev_kernel(t,0,&s,&ds,&c,&dc);
		r=sin_rev_round(2*c+(uint64_t)dc);
	}
	return reflected?pdi_neg(r):r;
}

void pdi_sin_cos_rev(pseudo_double_i x, pseudo_double_i* sin_out, pseudo_double_i* cos_out) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	bool negative,reflected;
	uint64_t t=sin_rev_fold(e,vx,&negative,&reflected);
	uint64_t s,c;
	int64_t ds,dc;
	sin_cos_rev_kernel(t,0,&s,&ds,&c,&dc);                              // one lookup serves both
	pseudo_double_i rs,rc;
	if(e<=-7) {
		rs=sin_rev_tiny(x,e,vx);
	} else {
		rs=t<PD_SIN_SMALL?sin_rev_small(t):sin_rev_round(2*s+(uint64_t)ds);
		rs=negative?pdi_neg(rs):rs;
	}
	rc=t>PD_QUARTER-PD_SIN_SMALL?sin_rev_small(PD_QUARTER-t):sin_rev_round(2*c+(uint64_t)dc);
	*sin_out=rs;
	*cos_out=reflected?pdi_neg(rc):rc;
}

pseudo_double_i pdi_tan_rev(pseudo_double_i x) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	bool negative,reflected;
	uint64_t t=sin_rev_fold(e,vx,&negative,&reflected);
	uint64_t s,c;
	int64_t ds,dc;
	sin_cos_rev_kernel(t,0,&s,&ds,&c,&dc);
	int32_t es=0,ec=0;
	signed_pd_internal vs,vc;
	if(e<=-7) {
		vs=x==0?0:sin_rev_tiny_u(e,vx<0?0-(uint64_t)vx:(uint64_t)vx,&es);
	} else {
		vs=t<PD_SIN_SMALL?sin_rev_small128_u(t,0,&es):sin_rev_round_u(2*s+(uint64_t)ds,&es);
	}
	vc=t>PD_QUARTER-PD_SIN_SMALL?sin_rev_small128_u(PD_QUARTER-t,0,&ec):sin_rev_round_u(2*c+(uint64_t)dc,&ec);
	return tan_rev_parts(vs,es,vc,ec,negative!=reflected);
}

// ---- radian versions: x is turned into turns with a 128-bit 1/(2*pi), so the reduction is accurate for |x| up to
// about 2^60; beyond that, too few fraction bits are left in x for the result to mean much anyway

// |x|/(2*pi) mod 1, as a 128-bit fraction of a turn (hi:lo): the caller handles the sign, so sin(-x) = -sin(x) exactly.
// 1/(2*pi) is used to 192 bits, which keeps the fraction exact enough for x close to a multiple of pi
static inline void sin_rad_to_turns(int32_t e, signed_pd_internal vx, uint64_t* hi, uint64_t* lo) {
	uint64_t ax=vx<0?0-(uint64_t)vx:(uint64_t)vx;
	// q2:q1:q0 = top 192 bits of |vx|*(1/(2*pi))*2^193
	uint64_t a1=multu64hi(ax,SIN_INV_TAU_0), a0=ax*SIN_INV_TAU_0;
	uint64_t b1=multu64hi(ax,SIN_INV_TAU_1), b0=ax*SIN_INV_TAU_1;
	uint64_t c1=multu64hi(ax,SIN_INV_TAU_2);
	uint64_t q0=c1+b0;
	uint64_t q1=a0+b1;
	uint64_t carry=(q1<b1);
	q1+=(q0<b0);
	carry+=(q1==0 && q0<b0);
	uint64_t q2=a1+carry;
	int32_t r=65-e;                                                    // |x|/(2*pi) = (q2:q1:q0)*2^(e-129): shift right by 65-e
	if(r<=-128 || r>=192) {
		*hi=*lo=0;                                                     // no fraction bits left, or a negligible fraction
	} else if(r<=-64) {
		*hi=q0<<(-r-64);
		*lo=0;
	} else if(r<0) {
		*hi=(q1<<-r)|(q0>>(64+r));
		*lo=q0<<-r;
	} else if(r==0) {
		*hi=q1;
		*lo=q0;
	} else if(r<64) {
		*hi=(q1>>r)|(q2<<(64-r));
		*lo=(q0>>r)|(q1<<(64-r));
	} else if(r==64) {
		*hi=q2;
		*lo=q1;
	} else if(r<128) {
		*hi=q2>>(r-64);
		*lo=(q1>>(r-64))|(q2<<(128-r));
	} else {
		*hi=0;
		*lo=q2>>(r-128);
	}
}

// as sin_rev_fold, on a 128-bit fraction of a turn for |x|; x_negative is the sign of x
static inline void sin_rad_fold(uint64_t* hi, uint64_t* lo, bool x_negative, bool* negative, bool* reflected) {
	*negative=(int64_t)*hi<0;                                          // second half turn
	if(*negative) {
		*hi=0-*hi-(*lo!=0);
		*lo=0-*lo;
	}
	*negative^=x_negative;                                             // sin(-x) = -sin(x); cos is even
	*reflected=*hi>(1ULL<<62) || (*hi==(1ULL<<62) && *lo!=0);
	if(*reflected) {
		*hi=(1ULL<<63)-*hi-(*lo!=0);
		*lo=0-*lo;
	}
}

// |sin(x radians)| for |x| <= 2^-7, before rounding; ax = |vx|. u = x*128/pi is computed directly, since the turn
// fraction would lose relative precision
static inline signed_pd_internal sin_rad_tiny_u(int32_t e, uint64_t ax, int32_t* exp_max) {
	uint64_t m=multu64hi(ax,(uint64_t)SIN_2_DIV_PI);                  // |u| = m*2^(e-57)
	uint64_t b=e<=-70?0:(m>>(-6-e));                                   // |u| * 2^63
	int64_t z=(int64_t)multu64hi(b,b);
	int32_t s=PD_CLZ(m)-1;
	*exp_max=PSEUDO_DOUBLE_EXP_BIAS+e+5-s;
	return mults64hi((int64_t)(m<<s),sin_poly(z,mults64hi(z,z)));
}

// sin(x radians) for |x| <= 2^-7: the magnitude is rounded, so sin(-x) = -sin(x) exactly
static inline pseudo_double_i sin_rad_tiny(pseudo_double_i x, int32_t e, signed_pd_internal vx) {
	if(x==0) {
		return 0;
	}
	int32_t exp_max;
	signed_pd_internal vr=sin_rad_tiny_u(e,vx<0?0-(uint64_t)vx:(uint64_t)vx,&exp_max);
	pseudo_double_i r=pseudo_round_helper(vr,exp_max,false);
	return vx<0?pdi_neg(r):r;
}

pseudo_double_i pdi_sin(pseudo_double_i x) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;           // x = vx*2^(e-64)
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	if(e<=-6) {
		return sin_rad_tiny(x,e,vx);
	}
	uint64_t t,tlo;
	bool negative,reflected;
	sin_rad_to_turns(e,vx,&t,&tlo);
	sin_rad_fold(&t,&tlo,vx<0,&negative,&reflected);
	pseudo_double_i r;
	if(t<(1ULL<<56)) {
		r=sin_rev_small128(t,tlo);
	} else {
		uint64_t s,c;
		int64_t ds,dc;
		sin_cos_rev_kernel(t,tlo,&s,&ds,&c,&dc);
		r=sin_rev_round(2*s+(uint64_t)ds);
	}
	return negative?pdi_neg(r):r;
}

pseudo_double_i pdi_cos(pseudo_double_i x) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	uint64_t t,tlo;
	bool negative,reflected;
	sin_rad_to_turns(e,vx,&t,&tlo);
	sin_rad_fold(&t,&tlo,vx<0,&negative,&reflected);
	uint64_t ct=(1ULL<<62)-t-(tlo!=0);                                 // 1/4-t, for cos(t) = sin(1/4-t)
	pseudo_double_i r;
	if(ct<(1ULL<<56)) {
		r=sin_rev_small128(ct,0-tlo);
	} else {
		uint64_t s,c;
		int64_t ds,dc;
		sin_cos_rev_kernel(t,tlo,&s,&ds,&c,&dc);
		r=sin_rev_round(2*c+(uint64_t)dc);
	}
	return reflected?pdi_neg(r):r;
}

void pdi_sin_cos(pseudo_double_i x, pseudo_double_i* sin_out, pseudo_double_i* cos_out) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	uint64_t t,tlo;
	bool negative,reflected;
	sin_rad_to_turns(e,vx,&t,&tlo);
	sin_rad_fold(&t,&tlo,vx<0,&negative,&reflected);
	uint64_t s,c;
	int64_t ds,dc;
	sin_cos_rev_kernel(t,tlo,&s,&ds,&c,&dc);                          // one lookup serves both
	pseudo_double_i r;
	if(e<=-6) {
		r=sin_rad_tiny(x,e,vx);
	} else {
		r=t<(1ULL<<56)?sin_rev_small128(t,tlo):sin_rev_round(2*s+(uint64_t)ds);
		r=negative?pdi_neg(r):r;
	}
	*sin_out=r;
	uint64_t ct=(1ULL<<62)-t-(tlo!=0);                                 // 1/4-t, for cos(t) = sin(1/4-t)
	r=ct<(1ULL<<56)?sin_rev_small128(ct,0-tlo):sin_rev_round(2*c+(uint64_t)dc);
	*cos_out=reflected?pdi_neg(r):r;
}

pseudo_double_i pdi_tan(pseudo_double_i x) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	uint64_t t,tlo;
	bool negative,reflected;
	sin_rad_to_turns(e,vx,&t,&tlo);
	sin_rad_fold(&t,&tlo,vx<0,&negative,&reflected);
	uint64_t s,c;
	int64_t ds,dc;
	sin_cos_rev_kernel(t,tlo,&s,&ds,&c,&dc);
	int32_t es=0,ec=0;
	signed_pd_internal vs,vc;
	if(e<=-6) {
		vs=x==0?0:sin_rad_tiny_u(e,vx<0?0-(uint64_t)vx:(uint64_t)vx,&es);
	} else {
		vs=t<(1ULL<<56)?sin_rev_small128_u(t,tlo,&es):sin_rev_round_u(2*s+(uint64_t)ds,&es);
	}
	uint64_t ct=(1ULL<<62)-t-(tlo!=0);                                 // 1/4-t, for cos(t) = sin(1/4-t)
	vc=ct<(1ULL<<56)?sin_rev_small128_u(ct,0-tlo,&ec):sin_rev_round_u(2*c+(uint64_t)dc,&ec);
	return tan_rev_parts(vs,es,vc,ec,negative!=reflected);
}

// BEGIN generated by gen_atan_table.py - do not edit
// entry k: atan(k/32) in eighth turns, as 1.63
static const uint64_t atan_table[33]={
	0ULL,366866921168144865ULL,733019165973179315ULL,1097750376323490839ULL,
	1460370590102087790ULL,1820213850616598198ULL,2176645144780590001ULL,2529066501404984965ULL,
	2876922122323524151ULL,3219702464204954152ULL,3556947234599589987ULL,3888247309078361832ULL,
	4213245614813964765ULL,4531637057867651645ULL,4843167595760423177ULL,5147632573313123190ULL,
	5444874448539492761ULL,5734780037346714272ULL,6017277402032338873ULL,6292332500319591573ULL,
	6559945700227971320ULL,6820148252640030787ULL,7072998799058243470ULL,7318579977620523033ULL,
	7556995176630566093ULL,7788365472136086064ULL,8012826774755461371ULL,8230527201165149541ULL,
	8441624677453327682ULL,8646284774869589253ULL,8844678773246323181ULL,9036981943387980342ULL,
	9223372036854775808ULL,
};
// atan(d) in eighth turns = u*G(w), u = 64*d, w = u^2; G = g0 + g1 w + ... + g4 w^4 (relative error ~2^-72.5)
#define ATAN_G0 (5871781006564002453LL)
#define ATAN_G1 (-1911387046407553LL)
#define ATAN_G2 (1119953347499LL)
#define ATAN_G3 (-781217393LL)
#define ATAN_G4 (593076LL)
#define ATAN_PI_DIV_4 (7244019458077122842LL) // pi/4 as 1.63: an eighth turn in radians
// END generated by gen_atan_table.py - do not edit

// u is 64*d as signed 1.63 (|d| < 2^-6); result is atan(d) in eighth turns, scaled by 2^67
static inline int64_t atan_poly(int64_t u) {
	int64_t w=mults64hi(u,u);                                          // u^2 * 2^62
	int64_t w2=mults64hi(w,w);                                         // u^4 * 2^60
	int64_t h=(ATAN_G2+mults64hi(ATAN_G3,w))+mults64hi(ATAN_G4,w2);    // 2^72
	return mults64hi(u,(ATAN_G0+mults64hi(ATAN_G1,w))+mults64hi(h,w2)); // u*G(w): G at 2^68
}

// atan(r) in eighth turns as 1.63: k = round(32*r), u = 64*d as signed 1.63, d = (r - k/32)/(1 + r*k/32)
static inline uint64_t atan_kernel(uint32_t k, int64_t u) {
	return atan_table[k]+(uint64_t)(atan_poly(u)>>4);
}

// x is a 2.62 unsigned fixed in the range [0,1]
// calculate atan_rev_64_fixed(x) = atan(x) in eighth turns
// result is 2.62 unsigned fixed in the range [0,1]
uint64_t atan_rev_64_fixed(uint64_t x) {
	uint32_t k=(uint32_t)(((x>>56)+1)>>1);                             // nearest 32nd
	int64_t n=(int64_t)(x-((uint64_t)k<<57));                          // x - k/32, as 2.62
	uint64_t dd=(1ULL<<62)+(x>>5)*k+(((x&31)*k)>>5);                    // 1 + x*k/32, as 2.62 (up to 2^63)
	return (atan_kernel(k,divs64hi((int64_t)((uint64_t)n<<4),(int64_t)(dd>>1)))+1)>>1;      // u*2^63 = 64*n/dd * 2^63; 1.63 -> 2.62 rounded
}

// add_const +- atan(a/b) in turns, for 0 <= a <= b, times k*2^(kexp-63) (k=0 for turns), with a single rounding.
// a = ma*2^(ea-BIAS-64) and b = mb*2^(eb-BIAS-64), with mantissas in [2^62,2^63) (ma = 0 for a = 0, except when
// add_const = 0). add_const is in eighth turns, scaled by 2^60; negative: atan(a/b) is subtracted. The mantissas may
// use all 63 bits (asin and acos pass a square root): a pdi mantissa's low bits are 0, which some shifts rely on
static inline pseudo_double_i atan2_core(uint64_t ma, int32_t ea, uint64_t mb, int32_t eb, int64_t add_const, int negative, int64_t k, int32_t kexp) {
	uint64_t t8=0;                                                     // atan(r) in eighth turns, 1.63
	uint32_t ki=0;
	int64_t q=0;                                                       // for r < 1/64: r = q*2^(ea-eb-62)
	if(ma!=0 && ea-eb>=-7) {
		uint32_t al=(uint32_t)(eb-ea);                                 // 0..7
		uint64_t ma2=ma>>al;                                           // a aligned to b's exponent (a floor if ma uses all its bits)
		uint32_t a26=(uint32_t)(ma2>>37), b26=(uint32_t)(mb>>37);     // top 26 bits: a26*32 + b26 still fits in 32 bits
		ki=(a26*32+(b26>>1)+2)/b26;                                    // round(32*r); the +2 keeps r >= 1/64 - tiny when ki >= 1
		if(ki!=0) {
			// 32*(a - b*k/32), exact: |n| < mb/2, so the bits above drop out. 32*a is shifted from ma itself, so a
			// full-length ma (from asin/acos) loses at most 2 bits here, rather than al of them
			int64_t n=(int64_t)((al<=5?ma<<(5-al):ma>>(al-5))-ki*mb);
			uint64_t dd=mb+(ma2>>5)*ki+(((ma2&31)*ki)>>5);             // b + a*k/32 (up to 2^64)
			t8=atan_kernel(ki,divs64hi(n>>1,(int64_t)(dd>>1)));       // u*2^63 = 64*d*2^63 = 2*n/dd * 2^63
		}
	}
	if(ki==0 && ma!=0) {
		// r < 1/64: r = a/b directly, and atan(r) = u*G(w) with u = 64*r
		q=divs64hi((int64_t)(ma>>2),(int64_t)mb);                      // q in (2^61,2^63)
		int32_t sh=ea-eb+7;                                            // u*2^63 = q*2^sh
		t8=atan_kernel(0,sh>0?(int64_t)((uint64_t)q<<sh):(sh>-64?q>>-sh:0));
	}
	signed_pd_internal vr;
	int32_t exp_max;
	if(add_const==0) {
		// q1: the result is atan(r), which can be small, so keep relative precision
		if(ki==0) {
			int32_t sh=ea-eb+7;
			int64_t u=sh>0?(int64_t)((uint64_t)q<<sh):(sh>-64?q>>-sh:0);
			int32_t s=PD_CLZ(q)-1;                                     // normalise r
			int64_t w=mults64hi(u,u);
			int64_t w2=mults64hi(w,w);
			int64_t h=(ATAN_G2+mults64hi(ATAN_G3,w))+mults64hi(ATAN_G4,w2);
			vr=mults64hi((int64_t)((uint64_t)q<<s),(ATAN_G0+mults64hi(ATAN_G1,w))+mults64hi(h,w2)); // turns = vr*2^(ea-eb-63-s)
			exp_max=PSEUDO_DOUBLE_EXP_BIAS+1+ea-eb-s;
		} else if(t8>>63) {
			vr=(signed_pd_internal)(t8>>1);                            // exactly 1 eighth turn
			exp_max=PSEUDO_DOUBLE_EXP_BIAS-1;
		} else {
			vr=(signed_pd_internal)t8;                                 // eighth turns * 2^63 = turns * 2^66
			exp_max=PSEUDO_DOUBLE_EXP_BIAS-2;
		}
	} else {
		// q2-q4: at least 1/8 turn, so absolute precision is enough
		int64_t t=(int64_t)(t8>>3);                                    // eighth turns * 2^60
		vr=add_const+(negative?-t:t);
		exp_max=PSEUDO_DOUBLE_EXP_BIAS+1;                              // eighth turns * 2^60 = turns * 2^63
	}
	if(k!=0) {
		int32_t z=PD_CLZ(vr)-1;                                        // normalise first, so the product keeps full precision
		vr=mults64hi((signed_pd_internal)((uint64_t)vr<<z),k);
		exp_max+=1+kexp-z;
	}
	return pseudo_round_helper(vr,exp_max,false);
}

// atan2(y,x) in turns, in (-1/2,1/2], times k*2^(kexp-63) (k=0 for turns), with a single rounding at the end
static inline pseudo_double_i pdi_atan2_scaled(pseudo_double_i y, pseudo_double_i x, int64_t k, int32_t kexp) {
	int negative=0;        // boolean: atan(r) is subtracted from add_const
	int64_t add_const;     // in eighth turns, scaled by 2^60
	int y_negative=0;      // boolean: q5-q8 are q4-q1 mirrored, so work with |y| and negate the result at the end
	if(((signed_pd_internal)y)<0) {
		y=pdi_neg(y);
		y_negative=1;
	}
	if(y==0) {
		if(((signed_pd_internal)x)>=0) {
			return 0;                                                  // includes atan2(0,0) = 0
		}
		add_const=4LL<<60;                                             // 1/2
	} else if(x==0) {
		add_const=2LL<<60;                                             // 1/4
		y=0;                                                           // so atan(r) = 0
	} else if(((signed_pd_internal)x)>0) {
		if(pdi_gte(x,y)) {
			// q1 (q8)
			add_const=0;
		} else {
			// q2 (q7)
			{ pseudo_double_i t=y; y=x; x=t; }
			add_const=2LL<<60;
			negative=1;
		}
	} else { // x<0
		x=pdi_neg(x);
		if(pdi_gte(x,y)) {
			// q4 (q5)
			add_const=4LL<<60;
			negative=1;
		} else {
			// q3 (q6)
			{ pseudo_double_i t=y; y=x; x=t; }
			add_const=2LL<<60;
		}
	}
	// now 0 <= y <= x, and the result is add_const +- atan(r) for r = y/x
	int32_t ea=(int32_t)(y&EXP_MASK), eb=(int32_t)(x&EXP_MASK);
	uint64_t ma=y==0?0:(uint64_t)(y&EXP_MASK_INV), mb=(uint64_t)(x&EXP_MASK_INV);
	pseudo_double_i r=atan2_core(ma,ea,mb,eb,add_const,negative,k,kexp); // rounds the magnitude, so atan2(-y,x) = -atan2(y,x)
	return y_negative?pdi_neg(r):r;
}

pseudo_double_i pdi_atan2_rev(pseudo_double_i y, pseudo_double_i x) {
	return pdi_atan2_scaled(y,x,0,0);
}

pseudo_double_i pdi_atan2(pseudo_double_i y, pseudo_double_i x) {
	return pdi_atan2_scaled(y,x,ATAN_PI_DIV_4,3);                      // turns * 2*pi = turns * 8 * pi/4
}

pseudo_double_i pdi_atan_rev(pseudo_double_i x) {
	return pdi_atan2_scaled(x,(1ULL<<62)+PSEUDO_DOUBLE_EXP_BIAS+2,0,0); // atan(x) = atan2(x,1)
}

pseudo_double_i pdi_atan(pseudo_double_i x) {
	return pdi_atan2_scaled(x,(1ULL<<62)+PSEUDO_DOUBLE_EXP_BIAS+2,ATAN_PI_DIV_4,3);
}

// ---- asin, acos ----

// sqrt(1-x^2) for |x| = ma*2^(ea-BIAS-64) (ma in [2^62,2^63)), as a full-length mantissa in [2^62,2^63) and exponent
// field, within a few units. 1-x^2 is formed exactly in 128 bits (except below 2^-126 for small x, where it
// doesn't matter), so x near 1 keeps full relative precision. Returns false for |x| > 1; *ms = 0 for |x| = 1
static inline bool pd_sqrt_1mx2(uint64_t ma, int32_t ea, uint64_t* ms, int32_t* es) {
	if(ea>=PSEUDO_DOUBLE_EXP_BIAS+2) {
		if(ea>PSEUDO_DOUBLE_EXP_BIAS+2 || ma>(1ULL<<62)) {
			return false;                                              // |x| > 1
		}
		*ms=0;                                                         // |x| = 1
		*es=0;
		return true;
	}
	// x^2 = (qh:ql)*2^(2*ea-2*BIAS-128); as 2.126 fixed that is a right shift by 2*(BIAS+1-ea) >= 0. Branch-free shift
	// by 0..127, as a shift by 0..63 and a 64-bit select (a floor, so T below is at most 2^-126 high)
	uint64_t qh=multu64hi(ma,ma), ql=ma*ma;
	int32_t d=2*(PSEUDO_DOUBLE_EXP_BIAS+1-ea);
	uint32_t n=d<127?(uint32_t)d:127;                                  // qh < 2^63, so 127 already clears it
	uint32_t s=n&63;
	uint64_t big=0-(uint64_t)(n>>6);
	uint64_t xl=(ql>>s)|((qh<<1)<<(63-s)), xh=qh>>s;                   // (qh<<1)<<(63-s) is qh<<(64-s), and 0 for s=0
	xl=(xl&~big)|(xh&big);
	xh&=~big;
	// T = 1-x^2 as 2.126, in (0,2^126]. For |x| < 1 it is at least about 2^(126-2*precision), so its high word is nonzero
	uint64_t tl=0-xl;
	uint64_t th=(1ULL<<62)-xh-(xl!=0);
	int32_t sh=PD_CLZ(th)&~1;                                          // normalise to [2^126,2^128) by an even shift
	th=(th<<sh)|((tl>>1)>>(63-sh));
	// sqrt(T) = v*2^63 with v in [1,2), from the top 64 bits: within [-4.1,+1.9] units of 2^-62. A Newton step using
	// all 128 bits would bring that to half a unit, but it costs about 4 ns and makes no measurable difference here
	uint64_t y=multu64hi(inv_sqrt64_fixed(th),th)<<1;                  // v as 2.62
	int32_t e=PSEUDO_DOUBLE_EXP_BIAS+2-(sh>>1);                        // sqrt(1-x^2) = y*2^(-62-sh/2)
	if(y>>63) {                                                        // v rounded up to 2
		y>>=1;
		e+=1;
	}
	*ms=y;
	*es=e;
	return true;
}

// asin(x) in turns, in [-1/4,1/4], times k*2^(kexp-63) (k=0 for turns): atan2(x, sqrt(1-x^2))
static inline pseudo_double_i pdi_asin_scaled(pseudo_double_i x, int64_t k, int32_t kexp) {
	if(x==0) {
		return 0;
	}
	int x_negative=((signed_pd_internal)x)<0;
	pseudo_double_i a=x_negative?pdi_neg(x):x;
	uint64_t ma=(uint64_t)(a&EXP_MASK_INV), ms;
	int32_t ea=(int32_t)(a&EXP_MASK), es;
	if(!pd_sqrt_1mx2(ma,ea,&ms,&es)) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_RANGE;                                             // |x| > 1
#endif
		return PD_NAN;
	}
	// |x| <= sqrt(1-x^2): atan(|x|/s), else 1/4 - atan(s/|x|) (which also gives exactly 1/4 for |x| = 1, as ms = 0)
	int q1=ms!=0 && (ea<es || (ea==es && ma<=ms));
	pseudo_double_i r=atan2_core(q1?ma:ms,q1?ea:es,q1?ms:ma,q1?es:ea,q1?0:(2LL<<60),!q1,k,kexp);
	return x_negative?pdi_neg(r):r;                                    // magnitude rounded, so asin(-x) = -asin(x)
}

// acos(x) in turns, in [0,1/2], times k*2^(kexp-63) (k=0 for turns): atan2(sqrt(1-x^2), x)
static inline pseudo_double_i pdi_acos_scaled(pseudo_double_i x, int64_t k, int32_t kexp) {
	if(x==0) {
		return atan2_core(0,0,1ULL<<62,PSEUDO_DOUBLE_EXP_BIAS+2,2LL<<60,0,k,kexp); // exactly 1/4
	}
	int x_negative=((signed_pd_internal)x)<0;
	pseudo_double_i a=x_negative?pdi_neg(x):x;
	uint64_t ma=(uint64_t)(a&EXP_MASK_INV), ms;
	int32_t ea=(int32_t)(a&EXP_MASK), es;
	if(!pd_sqrt_1mx2(ma,ea,&ms,&es)) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_RANGE;                                             // |x| > 1
#endif
		return PD_NAN;
	}
	if(ms==0 && !x_negative) {
		return 0;                                                      // acos(1)
	}
	// s <= |x|: x > 0 gives atan(s/|x|), x < 0 gives 1/2 - atan(s/|x|)
	// s > |x|:  x > 0 gives 1/4 - atan(|x|/s), x < 0 gives 1/4 + atan(|x|/s)
	int s_small=es<ea || (es==ea && ms<=ma);
	int64_t add_const=s_small?(x_negative?(4LL<<60):0):(2LL<<60);
	int negative=s_small?x_negative:!x_negative;
	return atan2_core(s_small?ms:ma,s_small?es:ea,s_small?ma:ms,s_small?ea:es,add_const,negative,k,kexp);
}

pseudo_double_i pdi_asin_rev(pseudo_double_i x) {
	return pdi_asin_scaled(x,0,0);
}

pseudo_double_i pdi_asin(pseudo_double_i x) {
	return pdi_asin_scaled(x,ATAN_PI_DIV_4,3);                         // turns * 2*pi = turns * 8 * pi/4
}

pseudo_double_i pdi_acos_rev(pseudo_double_i x) {
	return pdi_acos_scaled(x,0,0);
}

pseudo_double_i pdi_acos(pseudo_double_i x) {
	return pdi_acos_scaled(x,ATAN_PI_DIV_4,3);
}

// (+-a*b) + (+-c) with a single rounding (to nearest; ties to even when the sum is exact), all in 128 bits. The
// negations are done on the 128-bit forms, which have room for them, so the asymmetric two's complement range never
// gets in the way. pd_mul_add_core() only comes here for the cases its 64-bit fast path can't settle
static inline pseudo_double_i pd_mul_add_exact(pseudo_double_i a, pseudo_double_i b, pseudo_double_i c, bool neg_p, bool neg_c) {
	if(a==0 || b==0) {
		return neg_c?pdi_neg(c):c;
	}
	// p = va*vb exactly, as 128 bits: a*b = p*2^(ep-BIAS-128)
	signed_pd_internal va=(signed_pd_internal)(a&EXP_MASK_INV), vb=(signed_pd_internal)(b&EXP_MASK_INV);
	int64_t ph=mults64hi(va,vb);
	uint64_t pl=(uint64_t)va*(uint64_t)vb;
	if(neg_p) {                                                        // |p| <= 2^126, so -p fits
		pl=0-pl;
		ph=(int64_t)(~(uint64_t)ph+(pl==0));
	}
	int32_t ep=(int32_t)(a&EXP_MASK)+(int32_t)(b&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	int32_t e=ep;
	bool sticky=false;                                                 // bits were shifted out, so the sum is a floor
	if(c!=0) {
		// c as 128 bits, vc*2^62, so |p| <= 2^126 and |c| <= 2^125 add without overflow: c = (ch:cl)*2^(ec+2-BIAS-128)
		signed_pd_internal vc=(signed_pd_internal)(c&EXP_MASK_INV);
		int64_t ch=vc>>2;
		uint64_t cl=(uint64_t)vc<<62;
		if(neg_c) {                                                    // |c| <= 2^125 here, so -c fits
			cl=0-cl;
			ch=(int64_t)(~(uint64_t)ch+(cl==0));
		}
		int32_t ec=(int32_t)(c&EXP_MASK)+2;
		// shift the one with the smaller exponent right (a floor). Only |ep-ec| <= 1 can cancel heavily, and those
		// shifts are exact, because the low bits of p and c are 0
		int32_t d=ep-ec;
		int64_t* sh=d>=0?&ch:&ph;
		uint64_t* sl=d>=0?&cl:&pl;
		int32_t n=d>=0?d:-d;
		if(d<0) {
			e=ec;
		}
		if(n>=128) {
			sticky=(*sh!=0 || *sl!=0);
			*sl=(uint64_t)(*sh>>63);
			*sh=*sh>>63;
		} else if(n>=64) {
			sticky=(*sl!=0 || (n>64 && ((uint64_t)*sh<<(128-n))!=0));
			*sl=(uint64_t)(*sh>>(n-64));
			*sh=*sh>>63;
		} else if(n>0) {
			sticky=(*sl<<(64-n))!=0;
			*sl=(*sl>>n)|((uint64_t)*sh<<(64-n));
			*sh=*sh>>n;
		}
		uint64_t l=pl+cl;
		ph=(int64_t)((uint64_t)ph+(uint64_t)ch+(l<cl));
		pl=l;
	}
	if(ph==0 && pl==0) {
		return 0;                                                      // exact cancellation (a sticky sum can't be 0: see above)
	}
	// normalise the 128-bit sum (count its redundant sign bits), then round its top 64 bits once
	uint64_t sgn=(uint64_t)(ph>>63);
	uint64_t th=(uint64_t)ph^sgn, tl=pl^sgn;
	int32_t lz=th!=0?PD_CLZ(th)-1:(tl!=0?63+PD_CLZ(tl):127);
	uint64_t hi,lo;
	if(lz>=64) {
		hi=pl<<(lz-64);
		lo=0;
	} else if(lz>0) {
		hi=((uint64_t)ph<<lz)|(pl>>(64-lz));
		lo=pl<<lz;
	} else {
		hi=(uint64_t)ph;
		lo=pl;
	}
	sticky|=(lo!=0);
	return pseudo_round_helper((signed_pd_internal)hi,e-lz,!sticky);
}

// (+-a*b) + (+-c) with a single rounding. Fast path: the top 64 bits of a*b plus c in 64 bits, which is a floor
// that is less than 2 units low. It rounds the same way as the exact sum unless the guard bits are just under a
// midpoint, or the sum cancelled heavily; those (about 1.6% of random cases at 10 exponent bits) take the exact path
static inline pseudo_double_i pd_mul_add_core(pseudo_double_i a, pseudo_double_i b, pseudo_double_i c, bool neg_p, bool neg_c) {
	if(a==0 || b==0) {
		return neg_c?pdi_neg(c):c;
	}
	signed_pd_internal va=(signed_pd_internal)(a&EXP_MASK_INV), vb=(signed_pd_internal)(b&EXP_MASK_INV);
	int64_t p=mults64hi(va,vb);                                        // |p| <= 2^62, a floor
	uint64_t sticky=(uint64_t)va*(uint64_t)vb;                         // the low half
	if(neg_p) {
		p=-p-(sticky!=0);                                              // floor(-x) = -floor(x)-1 when inexact
	}
	int32_t ep=(int32_t)(a&EXP_MASK)+(int32_t)(b&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	int64_t vc=((signed_pd_internal)(c&EXP_MASK_INV))>>2;              // |c| <= 2^61, exact
	if(neg_c) {
		vc=-vc;
	}
	int32_t ec=(int32_t)(c&EXP_MASK)+2;                                // c = 0: exponent field 0, so it shifts away
	int32_t d=ep-ec;
	bool p_big=d>=0;
	int64_t bv=p_big?p:vc, sv=p_big?vc:p;
	int32_t e=p_big?ep:ec;
	uint32_t n=(uint32_t)(p_big?d:-d);
	n=n<63?n:63;
	sticky|=((uint64_t)sv<<1)<<(63-n);                                 // bits shifted out of the smaller one
	int64_t sum=bv+(sv>>n);                                            // < 2 units below the exact sum (if sticky)
	if(sticky==0) {
		return pseudo_round_helper(sum,e,true);                        // exact
	}
	if(sum!=0) {
		int32_t lb=PD_CLZ((uint64_t)(sum^(sum>>63))|1)-1;              // the helper's normalising shift (sum = -1: 62)
		uint64_t g=((uint64_t)sum<<lb)&((1ULL<<PSEUDO_DOUBLE_EXP_BITS)-1); // guard bits; the error is now < 2<<lb
		const uint64_t half=1ULL<<(PSEUDO_DOUBLE_EXP_BITS-1);
		if(lb<=PSEUDO_DOUBLE_EXP_BITS-2 && half-1-g>=(2ULL<<lb)-1) {   // [g,g+(2<<lb)) doesn't straddle the midpoint
			return pseudo_round_helper(sum,e,false);
		}
	}
	return pd_mul_add_exact(a,b,c,neg_p,neg_c);
}

// a*b + c with a single rounding
pseudo_double_i pdi_mul_add(pseudo_double_i a, pseudo_double_i b, pseudo_double_i c) {
	return pd_mul_add_core(a,b,c,false,false);
}

// a*b - c with a single rounding
pseudo_double_i pdi_mul_sub(pseudo_double_i a, pseudo_double_i b, pseudo_double_i c) {
	return pd_mul_add_core(a,b,c,false,true);
}

// c - a*b with a single rounding
pseudo_double_i pdi_sub_mul(pseudo_double_i a, pseudo_double_i b, pseudo_double_i c) {
	return pd_mul_add_core(a,b,c,true,false);
}

// ---- hypot, hypot3, inv_hypot, inv_hypot3 ----

// Adds x^2 to the unsigned 128-bit sum (*h:*l), with everything aligned to the exponent field emax (at least x's):
// the sum is then S*2^(2*emax-2*BIAS-128). The square is exact (up to 2^126); the alignment is a floor, and any bits
// shifted out are ORed into *st. Branch-free: the shift by 0..127 is split into a 64-bit select and a shift by 0..63
static inline void pd_add_square(pseudo_double_i x, int32_t emax, uint64_t* h, uint64_t* l, uint64_t* st) {
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	uint64_t ax=vx<0?0-(uint64_t)vx:(uint64_t)vx;                     // up to 2^63
	uint64_t qh=multu64hi(ax,ax), ql=ax*ax;
	int32_t d=2*(emax-(int32_t)(x&EXP_MASK));                          // x = 0 has exponent field 0: any shift works
	uint32_t dd=d<127?(uint32_t)d:127;                                 // qh < 2^63, so a shift of 127 already clears it
	uint32_t s=dd&63;
	uint64_t big=0-(uint64_t)(dd>>6);                                  // all ones for shifts of 64 or more
	uint64_t qh_out=(qh<<1)<<(63-s);                                   // qh<<(64-s), and 0 for s=0
	uint64_t lo=(ql>>s)|qh_out, hi=qh>>s;
	*st|=(((ql<<1)<<(63-s))&~big)|((ql|qh_out)&big);
	lo=(lo&~big)|(hi&big);
	hi&=~big;
	*l+=lo;
	*h+=hi+(*l<lo);
}

// Normalises a nonzero sum from pd_add_square to S in [2^126,2^128) by an even shift (exact: the largest square is at
// least 2^124, so the shift is 0 or 2). Returns half the shift
static inline int32_t pd_norm_squares(uint64_t* h, uint64_t* l) {
	int32_t sh=PD_CLZ(*h)&~1;
	*h=(*h<<sh)|((*l>>1)>>(63-sh));
	*l<<=sh;
	return sh>>1;
}

// sqrt(S*2^(2*emax-2*BIAS-128)) correctly rounded, for a sum from pd_add_square (sticky: bits were lost below S)
static inline pseudo_double_i pd_hypot_root(uint64_t s_hi, uint64_t s_lo, int32_t emax, bool sticky) {
	int32_t rexp=emax+1-pd_norm_squares(&s_hi,&s_lo);                 // result = v*2^63 * 2^(emax-BIAS-64-sh/2)
	// sqrt(S) = v*2^63 with v in [1,2): approximate v from the top 64 bits of S. y is within [-4.1,+1.9] units of v
	uint64_t y=multu64hi(inv_sqrt64_fixed(s_hi),s_hi)<<1;              // v as 2.62
	// q = nearest multiple of the result ulp, then checked exactly when y is close to a midpoint. Compare the midpoint m
	// with the root through m^2 (4.124) against v^2 = S/4: m^2 > v^2 when 4m^2 > S (with a sticky S the true value is
	// slightly above S, but 4m^2 and S are both multiples of the normalising shift, so the same test works); an exact
	// tie is possible here (hypot(3k,4k) = 5k)
	const uint64_t ulp=1ULL<<PSEUDO_DOUBLE_EXP_BITS, half=ulp>>1;
	uint64_t q=(y+half)&~(ulp-1);
	uint64_t frac=(y+half)&(ulp-1);
	if(frac<16 || frac>ulp-16) {
		for(int side=0;side<2;side++) {
			uint64_t m=side==0?q-half:q+half;
			uint64_t mh=multu64hi(m,m), ml=m*m;                        // m^2 (m < 2^63, so 4m^2 < 2^128)
			uint64_t fh=(mh<<2)|(ml>>62), fl=ml<<2;                    // 4m^2
			bool above=fh>s_hi || (fh==s_hi && fl>s_lo);               // m is above the root
			bool equal=fh==s_hi && fl==s_lo && !sticky;                // m is exactly the root: ties to even
			if(side==0) {
				if(above || (equal && (q&ulp))) {
					q-=ulp;                                            // lower midpoint is above the root (or a tie to an odd q)
					break;
				}
			} else if(!above && !(equal && !(q&ulp))) {
				q+=ulp;                                                // upper midpoint is below the root (or a tie to an odd q)
			}
		}
	}
	if(q==(1ULL<<63)) {                                                // rounded up to exactly 2
		q=1ULL<<62;
		rexp+=1;
	}
	return pseudo_round_helper((signed_pd_internal)q,rexp,false);      // q is already rounded: this only checks the range
}

// 1/sqrt(S*2^(2*emax-2*BIAS-128)), for a sum from pd_add_square
static inline pseudo_double_i pd_inv_hypot_root(uint64_t s_hi, uint64_t s_lo, int32_t emax) {
	int32_t rexp=2*PSEUDO_DOUBLE_EXP_BIAS+3-emax+pd_norm_squares(&s_hi,&s_lo);
	uint64_t y=inv_sqrt64_fixed(s_hi);                                 // 1/v as 1.63, in (1/2,1]
	return pseudo_round_helper((signed_pd_internal)(y>>1),rexp,false); // 1/sqrt(S) = y*2^-126
}

// sqrt(x^2+y^2), correctly rounded, with no overflow or underflow in between
pseudo_double_i pdi_hypot(pseudo_double_i x, pseudo_double_i y) {
	int32_t ex=(int32_t)(x&EXP_MASK), ey=(int32_t)(y&EXP_MASK);
	int32_t emax=ex>ey?ex:ey;
	uint64_t h=0,l=0,st=0;
	pd_add_square(x,emax,&h,&l,&st);
	pd_add_square(y,emax,&h,&l,&st);
	if(h==0) {
		return 0;                                                      // both 0 (otherwise the largest square is >= 2^124)
	}
	return pd_hypot_root(h,l,emax,st!=0);
}

// sqrt(x^2+y^2+z^2), correctly rounded
pseudo_double_i pdi_hypot3(pseudo_double_i x, pseudo_double_i y, pseudo_double_i z) {
	int32_t ex=(int32_t)(x&EXP_MASK), ey=(int32_t)(y&EXP_MASK), ez=(int32_t)(z&EXP_MASK);
	int32_t emax=ex>ey?ex:ey;
	emax=emax>ez?emax:ez;
	uint64_t h=0,l=0,st=0;
	pd_add_square(x,emax,&h,&l,&st);
	pd_add_square(y,emax,&h,&l,&st);
	pd_add_square(z,emax,&h,&l,&st);
	if(h==0) {
		return 0;
	}
	return pd_hypot_root(h,l,emax,st!=0);
}

// 1/sqrt(x^2+y^2)
pseudo_double_i pdi_inv_hypot(pseudo_double_i x, pseudo_double_i y) {
	int32_t ex=(int32_t)(x&EXP_MASK), ey=(int32_t)(y&EXP_MASK);
	int32_t emax=ex>ey?ex:ey;
	uint64_t h=0,l=0,st=0;
	pd_add_square(x,emax,&h,&l,&st);
	pd_add_square(y,emax,&h,&l,&st);
	if(h==0) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_RANGE;                                             // 1/0
#endif
		return PD_NAN;
	}
	return pd_inv_hypot_root(h,l,emax);
}

// 1/sqrt(x^2+y^2+z^2)
pseudo_double_i pdi_inv_hypot3(pseudo_double_i x, pseudo_double_i y, pseudo_double_i z) {
	int32_t ex=(int32_t)(x&EXP_MASK), ey=(int32_t)(y&EXP_MASK), ez=(int32_t)(z&EXP_MASK);
	int32_t emax=ex>ey?ex:ey;
	emax=emax>ez?emax:ez;
	uint64_t h=0,l=0,st=0;
	pd_add_square(x,emax,&h,&l,&st);
	pd_add_square(y,emax,&h,&l,&st);
	pd_add_square(z,emax,&h,&l,&st);
	if(h==0) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_RANGE;                                             // 1/0
#endif
		return PD_NAN;
	}
	return pd_inv_hypot_root(h,l,emax);
}

// BEGIN generated by gen_hyp_table.py - do not edit
// f = 1 + w*T(w), T(w) = t0 + t1 w + ... + t7 w^7, t_j stored at 2^(S+j)
#define SINH_S 67 // interpolation error in f ~2^-89.8
#define SINH_T0 (6148914691236517205LL)
#define SINH_T1 (153722867280912930LL)
#define SINH_T2 (1830034134296583LL)
#define SINH_T3 (12708570377060LL)
#define SINH_T4 (57766228986LL)
#define SINH_T5 (185148171LL)
#define SINH_T6 (440827LL)
#define SINH_T7 (813LL)
#define COSH_S 65 // interpolation error in f ~2^-85.5
#define COSH_T0 (4611686018427387904LL)
#define COSH_T1 (192153584101141163LL)
#define COSH_T2 (3202559735019019LL)
#define COSH_T3 (28594283348384LL)
#define COSH_T4 (158857129712LL)
#define COSH_T5 (601731556LL)
#define COSH_T6 (1653099LL)
#define COSH_T7 (3455LL)
#define ASINH_S 73 // interpolation error in f ~2^-93.7
#define ASINH_T0 (-6148914691236517205LL)
#define ASINH_T1 (21617278211378381LL)
#define ASINH_T2 (-100526777396663LL)
#define ASINH_T3 (534484819058LL)
#define ASINH_T4 (-3074806132LL)
#define ASINH_T5 (18632386LL)
#define ASINH_T6 (-117138LL)
#define ASINH_T7 (747LL)
// END generated by gen_hyp_table.py - do not edit

// T(w) for the small-argument polynomials (f = 1 + w*T(w)): w as 1.63 in [0,1), w2 = w^2*2^62, w4 = w^4*2^60, and
// t_j at 2^(S+j); result at 2^S. Estrin's scheme, laid out as log2_tail()
static inline int64_t hyp_poly(int64_t w, int64_t w2, int64_t w4, int64_t t0, int64_t t1, int64_t t2, int64_t t3, int64_t t4, int64_t t5, int64_t t6, int64_t t7) {
	int64_t a=t0+mults64hi(t1,w);                                      // 2^S
	int64_t b=t2+mults64hi(t3,w);                                      // 2^(S+2)
	int64_t c=t4+mults64hi(t5,w);                                      // 2^(S+4)
	int64_t d=t6+mults64hi(t7,w);                                      // 2^(S+6)
	return a+mults64hi(b,w2)+mults64hi(c+mults64hi(d,w2),w4);
}

// 1 + w*T(w) as 2.62, for T at 2^s
static inline uint64_t hyp_one_plus(int64_t w, int64_t t, int s) {
	return (1ULL<<62)+(uint64_t)(mults64hi(w,t)>>(s-63));
}

// e^|x| = s1*2^(ipart-63) for |x| >= 1/2 (|x| = va*2^(ex-64), va in [2^62,2^63)), from the exp2 kernel. Also gives
// e^-|x| = s2*2^(ipart-d-63) as s2>>d (0 when it is below 2^-64 of e^|x|). Returns false if e^|x| overflows
static inline bool hyp_exp_pair(uint64_t va, int32_t ex, uint64_t* s1, uint64_t* s2d, int32_t* ipart) {
	int64_t h=mults64hi((int64_t)va,EXP2_LOG2E_HI);                    // t = |x|*log2(e) in 128 bits, as in pdi_exp
	uint64_t l=va*(uint64_t)EXP2_LOG2E_HI;
	int64_t m=mulsu64hi((int64_t)va,EXP2_LOG2E_LO);
	pd_add128(&h,&l,m>>63,(uint64_t)m);
	uint64_t fraction=0;
	*ipart=0;
	if(exp2_split(h,l,ex+2,ipart,&fraction)==2) {                      // (|t| >= 0.72, so only 0 or 2 come back)
		return false;
	}
	*s1=exp2_64_fixed_inline(fraction);
	*s2d=0;
	if(*ipart<32) {
		// e^-t = 2^(-ipart-1) * 2^(1-fraction) (or 2^-ipart for a zero fraction): the same kernel on 1-fraction
		*s2d=exp2_64_fixed_inline(0-fraction)>>(2**ipart+(fraction!=0));
	}
	return true;
}

pseudo_double_i pdi_sinh(pseudo_double_i x) {
	if(x==0) {
		return 0;
	}
	int x_negative=((signed_pd_internal)x)<0;
	pseudo_double_i a=x_negative?pdi_neg(x):x;
	int32_t ex=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;           // |x| = va*2^(ex-64)
	uint64_t va=(uint64_t)(a&EXP_MASK_INV);
	signed_pd_internal vr;
	int32_t exp_max;
	if(ex<=0) {
		// |x| < 1/2: sinh(x) = x*(1 + w*T(w)), w = (2x)^2, which keeps full relative precision as x gets small
		int64_t u=ex<=-63?0:(int64_t)(va>>-ex);                        // 2|x| as 1.63
		int64_t w=mults64hi(u,u)<<1;
		int64_t w2=mults64hi(w,w), w4=mults64hi(w2,w2);
		uint64_t q=hyp_one_plus(w,hyp_poly(w,w2,w4,SINH_T0,SINH_T1,SINH_T2,SINH_T3,SINH_T4,SINH_T5,SINH_T6,SINH_T7),SINH_S);
		vr=(signed_pd_internal)multu64hi(va<<1,q);
		exp_max=ex+PSEUDO_DOUBLE_EXP_BIAS+1;
	} else {
		// sinh = (e^|x| - e^-|x|)/2: no worse than 1.6 bits of cancellation here
		uint64_t s1,s2d;
		int32_t ipart=0;
		if(!hyp_exp_pair(va,ex,&s1,&s2d,&ipart)) {
#if PD_ERROR_CHECK
			PD_DO_ERROR_OVERFLOW;
#endif
			return PD_NAN;
		}
		uint64_t dd=s1-s2d;                                            // sinh = dd*2^(ipart-64)
		uint32_t b=(uint32_t)(dd>>63);
		vr=(signed_pd_internal)(dd>>b);
		exp_max=PSEUDO_DOUBLE_EXP_BIAS+ipart+(int32_t)b;
	}
	pseudo_double_i r=pseudo_round_helper(vr,exp_max,false);
	return x_negative?pdi_neg(r):r;                                    // magnitude rounded, so sinh(-x) = -sinh(x)
}

pseudo_double_i pdi_cosh(pseudo_double_i x) {
	if(x==0) {
		return (1ULL<<62)+PSEUDO_DOUBLE_EXP_BIAS+2;
	}
	pseudo_double_i a=((signed_pd_internal)x)<0?pdi_neg(x):x;
	int32_t ex=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	uint64_t va=(uint64_t)(a&EXP_MASK_INV);
	if(ex<=0) {
		// |x| < 1/2: cosh(x) = 1 + w*T(w), w = (2x)^2
		int64_t u=ex<=-63?0:(int64_t)(va>>-ex);
		int64_t w=mults64hi(u,u)<<1;
		int64_t w2=mults64hi(w,w), w4=mults64hi(w2,w2);
		uint64_t q=hyp_one_plus(w,hyp_poly(w,w2,w4,COSH_T0,COSH_T1,COSH_T2,COSH_T3,COSH_T4,COSH_T5,COSH_T6,COSH_T7),COSH_S);
		return pseudo_round_helper((signed_pd_internal)q,PSEUDO_DOUBLE_EXP_BIAS+2,false);
	}
	// cosh = (e^|x| + e^-|x|)/2
	uint64_t s1,s2d;
	int32_t ipart=0;
	if(!hyp_exp_pair(va,ex,&s1,&s2d,&ipart)) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_OVERFLOW;
#endif
		return PD_NAN;
	}
	uint64_t ss=(s1>>1)+(s2d>>1);                                      // cosh = ss*2^(ipart-63)
	uint32_t b=(uint32_t)(ss>>63);
	return pseudo_round_helper((signed_pd_internal)(ss>>b),PSEUDO_DOUBLE_EXP_BIAS+ipart+1+(int32_t)b,false);
}

pseudo_double_i pdi_tanh(pseudo_double_i x) {
	if(x==0) {
		return 0;
	}
	int x_negative=((signed_pd_internal)x)<0;
	pseudo_double_i a=x_negative?pdi_neg(x):x;
	int32_t ex=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	uint64_t va=(uint64_t)(a&EXP_MASK_INV);
	signed_pd_internal vr;
	int32_t exp_max;
	if(ex<=0) {
		// |x| < 1/2: tanh(x) = x*(sinh(x)/x)/cosh(x) = x*(1 + d/cosh(x)), with d = sinh(x)/x - cosh(x) = w*(Ts - Tc)
		// formed from the two polynomial tails, so it keeps full precision
		int64_t u=ex<=-63?0:(int64_t)(va>>-ex);
		int64_t w=mults64hi(u,u)<<1;
		int64_t w2=mults64hi(w,w), w4=mults64hi(w2,w2);
		int64_t ts=hyp_poly(w,w2,w4,SINH_T0,SINH_T1,SINH_T2,SINH_T3,SINH_T4,SINH_T5,SINH_T6,SINH_T7);
		int64_t tc=hyp_poly(w,w2,w4,COSH_T0,COSH_T1,COSH_T2,COSH_T3,COSH_T4,COSH_T5,COSH_T6,COSH_T7);
		uint64_t qc=hyp_one_plus(w,tc,COSH_S);                         // cosh, 2.62
		int64_t d=mults64hi(w,(ts>>(SINH_S-COSH_S))-tc);               // at 2^(COSH_S-1), negative
		// |d/cosh| < 0.08 < 2^(62-COSH_S), so divs64hi's |d| < qc/2 holds: d/cosh at 2^(COSH_S+1)
		uint64_t q=(1ULL<<62)+(uint64_t)(divs64hi(d,(int64_t)qc)>>(COSH_S+1-62)); // tanh(x)/x, 2.62
		vr=(signed_pd_internal)multu64hi(va<<1,q);
		exp_max=ex+PSEUDO_DOUBLE_EXP_BIAS+1;
	} else {
		// tanh = 1 - 2/(e^2|x| + 1), from one exp2 kernel and one division
		int64_t h=mults64hi((int64_t)va,EXP2_LOG2E_HI);
		uint64_t l=va*(uint64_t)EXP2_LOG2E_HI;
		int64_t m=mulsu64hi((int64_t)va,EXP2_LOG2E_LO);
		pd_add128(&h,&l,m>>63,(uint64_t)m);
		int32_t ipart=0;
		uint64_t fraction=0;
		if(exp2_split(h,l,ex+3,&ipart,&fraction)==2 || ipart>=62) {   // t = 2|x|*log2(e) >= 1.44, so ipart >= 1
			vr=1LL<<62;                                                // 2/(e^2|x| + 1) < 2^-61: tanh rounds to 1
			exp_max=PSEUDO_DOUBLE_EXP_BIAS+2;
		} else {
			uint64_t s=exp2_64_fixed_inline(fraction);                 // e^2|x| = s*2^(ipart-63)
			int64_t dm=(int64_t)((s>>2)+((1ULL<<61)>>ipart));          // (e^2|x| + 1)*2^(61-ipart), < 2^62+2^60
			uint64_t q=(uint64_t)divs64hi((int64_t)(1ULL<<(61-ipart)),dm); // 2/(e^2|x| + 1)*2^63, at least 3
			vr=(signed_pd_internal)((1ULL<<63)-q);                     // tanh*2^63
			exp_max=PSEUDO_DOUBLE_EXP_BIAS+1;
		}
	}
	pseudo_double_i r=pseudo_round_helper(vr,exp_max,false);
	return x_negative?pdi_neg(r):r;
}

// ---- inverse hyperbolics ----

// sqrt of the unsigned 128-bit H = (h:l) > 0: returns y in [2^62,2^63) with sqrt(H) = y*2^*k, to about half a unit
static inline uint64_t hyp_sqrt128(uint64_t h, uint64_t l, int32_t* k) {
	int32_t sh;                                                        // normalise to [2^126,2^128) by an even shift
	if(h!=0) {
		sh=PD_CLZ(h)&~1;
		h=(h<<sh)|((l>>1)>>(63-sh));
		l<<=sh;
	} else {
		sh=(64+PD_CLZ(l))&~1;
		h=l<<(sh-64);
		l=0;
	}
	// sqrt(H) = v*2^63, v in [1,2): from the top 64 bits (within [-4.1,+1.9] units of 2^-62), then one Newton step
	// using all 128 bits: v += (H/4 - v^2)/(2v), with 1/v from the same inverse square root
	uint64_t w=inv_sqrt64_fixed(h);                                    // 1/v as 1.63
	uint64_t y=multu64hi(w,h)<<1;                                      // v as 2.62
	uint64_t yh=multu64hi(y,y), yl=y*y;
	uint64_t t4h=h>>2, t4l=(l>>2)|(h<<62);                             // H/4, on the same scale as y^2
	uint64_t rl=t4l-yl;
	uint64_t rh=t4h-yh-(t4l<yl);                                       // R = H/4 - y^2, |R| < 2^67
	int64_t r16=(int64_t)((rh<<60)|(rl>>4));                           // R/16 fits in 64 bits
	y+=(uint64_t)((mulsu64hi(r16,w)+(1LL<<57))>>58);                   // R/(2y) = R*w*2^-126 = R/16*w*2^-122, rounded
	*k=1-(sh>>1);                                                      // sqrt(H) = 2y*2^(-sh/2)
	if(y>>63) {                                                        // v rounded up to 2
		y>>=1;
		*k+=1;
	}
	return y;
}

// v*2^n as a 128-bit unsigned (hi:lo), for n in [-63,63] (a floor for n < 0)
static inline void hyp_u128(uint64_t v, int32_t n, uint64_t* hi, uint64_t* lo) {
	if(n>=0) {
		*hi=n==0?0:v>>(64-n);
		*lo=v<<n;
	} else {
		*hi=0;
		*lo=v>>-n;
	}
}

// log2 of the 64.64 fixed value z = (hi:lo) >= 1, before rounding (as log2_general_u)
static inline signed_pd_internal hyp_log2_z(uint64_t hi, uint64_t lo, int32_t* exp_max) {
	int32_t c=PD_CLZ(hi);                                              // z = 2^e*(1+f)
	uint64_t f=c==63?lo:(hi<<(c+1))|(lo>>(63-c));
	return log2_general_u(63-c,f,exp_max);
}

// log2 result (vr, exp_max) times ln(2), rounded once
static inline pseudo_double_i hyp_ln_round(signed_pd_internal vr, int32_t exp_max) {
	int32_t z=PD_CLZ(vr>0?vr:~vr)-1;
	vr=mults64hi((signed_pd_internal)((uint64_t)vr<<z),LOG2_TO_LN);
	return pseudo_round_helper(vr,exp_max+1-z,false);
}

// asinh(x) = ln(|x| + sqrt(x^2+1)), with the sign of x
pseudo_double_i pdi_asinh(pseudo_double_i x) {
	if(x==0) {
		return 0;
	}
	int x_negative=((signed_pd_internal)x)<0;
	pseudo_double_i a=x_negative?pdi_neg(x):x;
	int32_t ex=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;           // |x| = va*2^(ex-64), in [2^(ex-2),2^(ex-1))
	uint64_t va=(uint64_t)(a&EXP_MASK_INV);
	pseudo_double_i r;
	if(ex<=-3) {
		// |x| < 1/16: asinh(x) = x*(1 + w*T(w)), w = (16x)^2
		int64_t u=ex<=-66?0:(int64_t)(va>>(-3-ex));                    // 16|x| as 1.63
		int64_t w=mults64hi(u,u)<<1;
		int64_t w2=mults64hi(w,w), w4=mults64hi(w2,w2);
		uint64_t q=hyp_one_plus(w,hyp_poly(w,w2,w4,ASINH_T0,ASINH_T1,ASINH_T2,ASINH_T3,ASINH_T4,ASINH_T5,ASINH_T6,ASINH_T7),ASINH_S);
		r=pseudo_round_helper((signed_pd_internal)multu64hi(va<<1,q),ex+PSEUDO_DOUBLE_EXP_BIAS+1,false);
	} else {
		signed_pd_internal vr;
		int32_t exp_max;
		if(ex>=34) {
			// |x| >= 2^32: ln(2|x|), as the rest is below 2^-66
			vr=log2_general_u(ex-1,va<<2,&exp_max);                    // |x| = 2^(ex-2)*(1+f)
		} else {
			// H = x^2+1 exactly (to 2^-126 for |x| < 1), z = |x| + sqrt(H) as 64.64, then log2(z)
			uint64_t qh=multu64hi(va,va), ql=va*va;                    // x^2 = (qh:ql)*2^-dd
			int32_t dd=2*(64-ex);                                      // 62..132
			int32_t scale;
			if(dd<=126) {
				if(dd>=64) {
					qh+=1ULL<<(dd-64);
				} else {
					ql+=1ULL<<dd;
					qh+=ql<(1ULL<<dd);                                 // carry
				}
				scale=-dd;
			} else {
				int32_t s=dd-126;                                      // 2, 4 or 6
				ql=(ql>>s)|(qh<<(64-s));
				qh=(qh>>s)+(1ULL<<62);
				scale=-126;
			}
			int32_t k;
			uint64_t y=hyp_sqrt128(qh,ql,&k);                          // sqrt(x^2+1) = y*2^(k+scale/2), >= 1
			uint64_t zh,zl,sh,sl;
			hyp_u128(va,ex,&zh,&zl);                                   // |x|*2^64 (exact)
			hyp_u128(y,k+scale/2+64,&sh,&sl);
			zl+=sl;
			zh+=sh+(zl<sl);
			vr=hyp_log2_z(zh,zl,&exp_max);
		}
		r=hyp_ln_round(vr,exp_max);
	}
	return x_negative?pdi_neg(r):r;
}

// acosh(x) = ln(x + sqrt(x^2-1)) for x >= 1
pseudo_double_i pdi_acosh(pseudo_double_i x) {
	int32_t ex=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	if(vx<=0 || ex<2) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_RANGE;                                             // x < 1
#endif
		return PD_NAN;
	}
	uint64_t m=(uint64_t)vx;
	if(ex==2 && m==(1ULL<<62)) {
		return 0;                                                      // acosh(1)
	}
	uint64_t tm=m-(1ULL<<62);                                          // for x < 2: x-1 = tm*2^-62 (exact)
	if(ex==2 && tm<(1ULL<<55)) {
		// x-1 < 2^-7: acosh(x) = 2*asinh(z), z = sqrt((x-1)/2) < 1/16, from the asinh polynomial (full relative
		// precision as x approaches 1, where the log form would need it from a value near 1)
		int32_t k;
		uint64_t y=hyp_sqrt128(0,tm<<1,&k);                            // sqrt(2*tm) = y*2^k, so z = y*2^(k-32)
		int32_t kz=k-32;                                               // z < 1/16, so kz <= -67
		int64_t u=(int64_t)(y>>(-67-kz));                              // 16z as 1.63
		int64_t w=mults64hi(u,u)<<1;
		int64_t w2=mults64hi(w,w), w4=mults64hi(w2,w2);
		uint64_t q=hyp_one_plus(w,hyp_poly(w,w2,w4,ASINH_T0,ASINH_T1,ASINH_T2,ASINH_T3,ASINH_T4,ASINH_T5,ASINH_T6,ASINH_T7),ASINH_S);
		return pseudo_round_helper((signed_pd_internal)multu64hi(y<<1,q),kz+66+PSEUDO_DOUBLE_EXP_BIAS,false); // 2*z*q*2^-62
	}
	signed_pd_internal vr;
	int32_t exp_max;
	if(ex>=34) {
		vr=log2_general_u(ex-1,m<<2,&exp_max);                         // x >= 2^32: ln(2x), as the rest is below 2^-66
	} else {
		// H = x^2-1 exactly: x^2 = m^2*2^-dd with dd = 2*(64-ex) in [62,124], and 1 = 2^dd on that scale
		int32_t dd=2*(64-ex);
		uint64_t qh=multu64hi(m,m), ql=m*m;
		if(dd>=64) {
			qh-=1ULL<<(dd-64);
		} else {
			qh-=(ql<(1ULL<<dd));
			ql-=1ULL<<dd;
		}
		int32_t k;
		uint64_t y=hyp_sqrt128(qh,ql,&k);                              // sqrt(x^2-1) = y*2^(k-dd/2), at least 1/8
		uint64_t zh,zl,sh,sl;
		hyp_u128(m,ex,&zh,&zl);                                        // x*2^64 (exact)
		hyp_u128(y,k-dd/2+64,&sh,&sl);
		zl+=sl;
		zh+=sh+(zl<sl);
		vr=hyp_log2_z(zh,zl,&exp_max);                                 // log2(x + sqrt(x^2-1)), at least 0.18
	}
	return hyp_ln_round(vr,exp_max);
}

// atanh(x) = (log2(1+x) - log2(1-x))*ln(2)/2 for |x| < 1: the two logs have opposite signs, so nothing cancels
pseudo_double_i pdi_atanh(pseudo_double_i x) {
	if(x==0) {
		return 0;
	}
	int x_negative=((signed_pd_internal)x)<0;
	pseudo_double_i a=x_negative?pdi_neg(x):x;
	int32_t ex=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	uint64_t va=(uint64_t)(a&EXP_MASK_INV);
	if(ex>=2) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_RANGE;                                             // |x| >= 1
#endif
		return PD_NAN;
	}
	signed_pd_internal vr;
	int32_t exp_max;
	if(ex<=-5) {
		// |x| < 2^-6: both logs from log2_small_u (relative precision), on the same scale
		int32_t s=-5-ex;
		int64_t rr=ex<=-68?0:(int64_t)(va>>s);
		vr=log2_small_u((int64_t)va,s,rr,&exp_max)-log2_small_u(-(int64_t)va,s,-rr,&exp_max);
	} else {
		uint64_t xf=ex>=0?va<<ex:va>>-ex;                              // |x| as 0.64 (exact)
		int32_t em1,em2;
		signed_pd_internal v1=log2_general_u(0,xf,&em1);               // log2(1+|x|)
		uint64_t om=0-xf;                                              // 1-|x| as 0.64 (exact, nonzero)
		int32_t c=PD_CLZ(om);
		signed_pd_internal v2=log2_general_u(-1-c,c==63?0:om<<(c+1),&em2); // log2(1-|x|) < 0
		exp_max=(em1>em2?em1:em2)+1;                                   // a bit of headroom for the sum
		vr=(v1>>(exp_max-em1))-(v2>>(exp_max-em2));
	}
	pseudo_double_i r=hyp_ln_round(vr,exp_max-1);                      // times ln(2)/2
	return x_negative?pdi_neg(r):r;
}

// BEGIN generated by gen_cbrt_table.py - do not edit
#define CBRT_C1_SHIFT 44 // c1 is stored as c1*2^20
#define CBRT_C2_SHIFT 40 // c2 is stored as c2*2^24
// initial approximation error ~2^-16.7
static const uint64_t inv_cbrt_packed[56]={ // c0:32 | c1:16 | c2:16
	0xFFFF67A4A9FEC5DFULL,0xF6246FD391679897ULL,0xEDA5F1247E6B78C8ULL,0xE63790B06F6161AAULL,
	0xDFA2CEC063365063ULL,0xD9BF88DD5931432AULL,0xD46F69E050D138D8ULL,0xCF9B037B49BA30A6ULL,
	0xCB2FE4E643A82A0BULL,0xC71F4E573E6924A6ULL,0xC35D48FE39D62031ULL,0xBFE0015E35D21C78ULL,
	0xBC9F4E9832441955ULL,0xB99458E82F1A16ABULL,0xB6B956302C451463ULL,0xB409566229BA126BULL,
	0xB1801B8C276D10B6ULL,0xAF19FA8D25570F39ULL,0xACD3C24323700DEAULL,0xAAAAA7BC21B30CC3ULL,
	0xA89C3630201B0BBEULL,0xA6A641F91EA30AD5ULL,0xA4C6DDE51D490A06ULL,0xA2FC52601C08094CULL,
	0xA145161F1ADF08A4ULL,0x9F9FC7FC19CB080DULL,0x9E0B29C918C90785ULL,0x9C861BF117D90708ULL,
	0x9B0F99C216F80697ULL,0x99A6B63E16250630ULL,0x984A9961155F05D2ULL,0x96FA7DC714A5057BULL,
	0x95B5AEA213F5052CULL,0x947B85F8135004E3ULL,0x934B6B1612B3049FULL,0x9224D13712200461ULL,
	0x9107365511940427ULL,0x8FF22221110F03F1ULL,0x8EE52510109103BFULL,0x8DDFD79310190391ULL,
	0x8CE1D9570FA70365ULL,0x8BEAD0A20F3A033DULL,0x8AFA69C00ED30317ULL,0x8A1056800E7002F4ULL,
	0x892C4DBB0E1102D3ULL,0x884E0AEC0DB702B3ULL,0x87754DD00D610296ULL,0x86A1DA0E0D0E027BULL,
	0x85D376EC0CBF0261ULL,0x8509EF050C720248ULL,0x8445100A0C290231ULL,0x8384AA880BE3021CULL,
	0x82C891B40BA00207ULL,0x82109B370B5F01F3ULL,0x815C9F050B2101E1ULL,0x80AC77340AE401D0ULL,
};
// END generated by gen_cbrt_table.py - do not edit

// x^(-1/3) for x as 3.61 in [1,8), as 0.64 in (0.5,1) (x = 1 gives 2^64-1)
static inline uint64_t inv_cbrt64_fixed_inline(uint64_t x) {
	// initial approximation: a quadratic on one of 56 segments of [1,8), good to about 2^-16.7
	uint64_t w=inv_cbrt_packed[(x>>58)-8];                        // segment from the top 6 bits; c0:32 | c1:16 | c2:16
	uint64_t u=x<<6;                                              // position within the segment, 0.64
	uint64_t c1=((w>>16)&0xFFFF)<<CBRT_C1_SHIFT;                  // 0.64
	uint64_t c2=(w&0xFFFF)<<CBRT_C2_SHIFT;                        // 0.64
	uint64_t y=(w&0xFFFFFFFF00000000ULL)-multu64hi(u,c1-multu64hi(u,c2)); // 0.64
	// two Newton-Raphson steps, y=y*(4-x*y^3)/3: error ~2^-16.7 -> ~2^-32 -> ~2^-63. The error term is -2e^2, so each
	// step lands below the true value, and the first one stays under 1 (so fits 0.64). (4-x*y^3)*2^61 is ~3*2^61, and
	// times 2/3 (then doubled) gives (4-x*y^3)/3 as 1.63
	uint64_t k=multu64hi((1ULL<<63)-multu64hi(x,multu64hi(multu64hi(y,y),y)),0xAAAAAAAAAAAAAAABULL)<<1;
	y=multu64hi(y,k)<<1;
	k=multu64hi((1ULL<<63)-multu64hi(x,multu64hi(multu64hi(y,y),y)),0xAAAAAAAAAAAAAAABULL)<<1;
	y=multu64hi(y,k);                                             // 1.63
	return y>>63?~0ULL:y<<1;                                      // 0.64: the result is 1 only for x = 1
}

// x is a 3.61 unsigned fixed in the range [1,8)
// result is x^(-1/3) as 1.63 unsigned fixed in the range (0.5,1]
uint64_t inv_cbrt64_fixed(uint64_t x) {
	return (inv_cbrt64_fixed_inline(x)>>1)+(x==(1ULL<<61)?1:0);   // exactly 1 for x = 1
}

// compares mid^3 with t = (t2:t1:0) (192 bits): >0 if mid^3 is bigger, <0 if smaller (never equal here)
static inline int cbrt_cube_cmp(uint64_t mid, uint64_t t2, uint64_t t1) {
	uint64_t h=multu64hi(mid,mid), l=mid*mid;                     // mid^2
	uint64_t c0=mid*l;
	uint64_t c1a=multu64hi(mid,l), c1=c1a+mid*h;
	uint64_t c2=multu64hi(mid,h)+(c1<c1a);                        // mid^3 = (c2:c1:c0)
	if(c2!=t2) {
		return c2>t2?1:-1;
	}
	if(c1!=t1) {
		return c1>t1?1:-1;
	}
	return c0!=0;
}

// cube root, correctly rounded (cbrt(-x) = -cbrt(x))
pseudo_double_i pdi_cbrt(pseudo_double_i x) {
	if(x==0) {
		return 0;
	}
	int x_negative=((signed_pd_internal)x)<0;
	pseudo_double_i a=x_negative?pdi_neg(x):x;
	int32_t e=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS-2;          // |x| = (ax/2^62)*2^e, ax in [2^62,2^63)
	uint64_t ax=(uint64_t)(a&EXP_MASK_INV);
	// e = 3q + r with r in {0,1,2}: |x| = m*2^(3q), m = ax*2^(r-62) in [1,8), cbrt(|x|) = cbrt(m)*2^q
	int32_t q=(e>=0?e:e-2)/3;                                          // floor(e/3)
	int32_t r=e-3*q;
	uint64_t y=inv_cbrt64_fixed_inline(r==0?ax>>1:ax<<(r-1));          // m^(-1/3) as 0.64, with m as 3.61 (exact: low bits of ax are 0)
	uint64_t t=multu64hi(ax,y)<<r;                                     // m*y*2^(62-r) * 2^r = m^(2/3)*2^62, in [2^62,2^64)
	uint64_t v=multu64hi(t,y);                                         // m*y^2*2^62: cbrt(m) as 2.62, within a few units
	// round to the nearest multiple of the ulp, then check that exactly when v is close to a midpoint: compare the
	// midpoint with the root through cubes, mid^3 against m*2^186 = ax*2^(124+r). An exact tie is impossible (mid^3 has
	// 3*(EXP_BITS-1) trailing zero bits, ax*2^(124+r) at least EXP_BITS+124)
	const uint64_t ulp=1ULL<<PSEUDO_DOUBLE_EXP_BITS, half=ulp>>1;
	uint64_t qr=(v+half)&~(ulp-1);
	uint64_t frac=(v+half)&(ulp-1);
	if(frac<32 || frac>ulp-32) {                                       // v is within [-8.8,+6.3] units of the root
		uint64_t t2=ax>>(4-r), t1=ax<<(60+r);
		if(cbrt_cube_cmp(qr-half,t2,t1)>0) {
			qr-=ulp;                                                   // lower midpoint is above the root
		} else if(cbrt_cube_cmp(qr+half,t2,t1)<0) {
			qr+=ulp;                                                   // upper midpoint is below the root
		}
	}
	if(qr==(1ULL<<63)) {                                               // rounded up to exactly 2
		qr=1ULL<<62;
		q+=1;
	}
	pseudo_double_i res=qr+(uint64_t)(q+2+PSEUDO_DOUBLE_EXP_BIAS);     // the exponent always stays in range
	return x_negative?pdi_neg(res):res;
}

// ---- powi: x^n for integer n ----
// Binary exponentiation on 128-bit mantissas (h:l) normalised to [2^127,2^128), value (h:l)*2^e. Every product is
// formed exactly and truncated to 128 bits (a floor), with a sticky flag for anything dropped. Each step loses under
// 2^-127, so even |n| = 2^31 (62 steps) stays within 2^-121 of the true value: the result rounds correctly unless
// the true value lies within that of a rounding midpoint. An exact result (nothing dropped) rounds ties to even

// Normalises a product (p3:p2:p1:p0) in [2^254,2^256) to its top 128 bits in [2^127,2^128) (branch-free), ORs the
// dropped bits into *sticky, and returns the exponent to add (128, or 127 if it took one more bit)
static inline int32_t powi_norm(uint64_t p3, uint64_t p2, uint64_t p1, uint64_t p0, uint64_t* h, uint64_t* l, uint64_t* sticky) {
	uint32_t s=(uint32_t)(1-(p3>>63));                                 // 1 if the product is in [2^254,2^255)
	*h=(p3<<s)|((p2>>1)>>(63-s));
	*l=(p2<<s)|((p1>>1)>>(63-s));
	*sticky|=(p1<<s)|p0;
	return 128-(int32_t)s;
}

// (h:l)^2 exactly: h^2*2^128 + 2*h*l*2^64 + l^2
static inline int32_t powi_sqr(uint64_t* h, uint64_t* l, uint64_t* sticky) {
	uint64_t ah=*h, al=*l;
	uint64_t p3=multu64hi(ah,ah), p2=ah*ah;
	uint64_t mh=multu64hi(ah,al), ml=ah*al;                            // h*l, doubled below (it is under 2^128)
	uint64_t dh=(mh<<1)|(ml>>63), dl=ml<<1;
	p3+=mh>>63;
	uint64_t p0=al*al, p1=multu64hi(al,al)+dl;
	uint64_t c=p1<dl;
	p2+=dh;
	p3+=p2<dh;
	p2+=c;
	p3+=p2<c;
	return powi_norm(p3,p2,p1,p0,h,l,sticky);
}

// (h:l)*(bh:bl) exactly
static inline int32_t powi_mul(uint64_t* h, uint64_t* l, uint64_t bh, uint64_t bl, uint64_t* sticky) {
	uint64_t ah=*h, al=*l;
	uint64_t p3=multu64hi(ah,bh), p2=ah*bh;
	uint64_t x1h=multu64hi(ah,bl), x1l=ah*bl;
	uint64_t x2h=multu64hi(al,bh), x2l=al*bh;
	uint64_t p0=al*bl;
	uint64_t p1=multu64hi(al,bl)+x1l;
	uint64_t c=p1<x1l;                                                 // carries out of word 1
	p1+=x2l;
	c+=p1<x2l;
	p2+=x1h;
	p3+=p2<x1h;
	p2+=x2h;
	p3+=p2<x2h;
	p2+=c;
	p3+=p2<c;
	return powi_norm(p3,p2,p1,p0,h,l,sticky);
}

// (h:l)*bh exactly, for a base whose low word is 0 (|x| itself): a 192-bit product
static inline int32_t powi_mul64(uint64_t* h, uint64_t* l, uint64_t bh, uint64_t* sticky) {
	uint64_t ah=*h, al=*l;
	uint64_t p3=multu64hi(ah,bh), p2=ah*bh;
	uint64_t p1=al*bh, x=multu64hi(al,bh);
	p2+=x;
	p3+=p2<x;
	return powi_norm(p3,p2,p1,0,h,l,sticky);
}

// floor(r*2^63/d) for r < d, d in (2^62,2^63), with the remainder: divs64hi needs its |x| < |y|/2, so divide
// floor(r/2)*2^64 and add back the last bit of r
static inline uint64_t powi_div63(uint64_t r, uint64_t d, uint64_t* rem) {
	uint64_t q=(uint64_t)divs64hi((int64_t)(r>>1),(int64_t)d);
	uint64_t s=(0-q*d)+((r&1)<<63);                                    // remainder so far (< d) + the last bit: < 3d
	uint32_t extra=(uint32_t)(s>=d)+(uint32_t)(s>=2*d);
	*rem=s-extra*d;
	return q+extra;
}

pseudo_double_i pdi_powi(pseudo_double_i x, int32_t n) {
	if(n==0) {
		return (1ULL<<62)+PSEUDO_DOUBLE_EXP_BIAS+2;                    // including 0^0 = 1
	}
	if(n==1) {
		return x;
	}
	if(n==2) {
		// x^2 = vx^2*2^(2*ex-128): the 128-bit square is exact, so the low word says whether anything is lost
		signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
		int32_t ex=(int32_t)(x&EXP_MASK);
		return pseudo_round_helper(mults64hi(vx,vx),2*ex-PSEUDO_DOUBLE_EXP_BIAS,(uint64_t)vx*(uint64_t)vx==0);
	}
	if(x==0) {
		if(n>0) {
			return 0;
		}
#if PD_ERROR_CHECK
		PD_DO_ERROR_OVERFLOW;                                          // 0 to a negative power
#endif
		return PD_NAN;
	}
	uint32_t un=n<0?0u-(uint32_t)n:(uint32_t)n;                        // |n|, also for n = INT32_MIN
	int negative=((signed_pd_internal)x)<0 && (un&1);
	pseudo_double_i a=((signed_pd_internal)x)<0?pdi_neg(x):x;
	int32_t ex=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;           // |x| = ax*2^(ex-64), ax in [2^62,2^63)
	uint64_t ax=(uint64_t)(a&EXP_MASK_INV);
	// the base b = |x| or 1/|x|, as (bh:bl)*2^be
	uint64_t bh,bl,sticky=0;
	int64_t be;
	if(n>0) {
		bh=ax<<1;                                                      // |x| = (ax*2^65)*2^(ex-129)
		bl=0;
		be=ex-129;
	} else if(ax==(1ULL<<62)) {
		bh=1ULL<<63;                                                   // |x| = 2^(ex-2): 1/|x| = 2^127*2^(-125-ex), exactly
		bl=0;
		be=-125-ex;
	} else {
		// 1/|x| = 2^(64-ex)/ax, from R = floor(2^188/ax) by long division: 2^125/ax = q1 + rem1/ax, then
		// R = q1*2^63 + floor(rem1*2^63/ax). ax > 2^62, so q1 < 2^63 and divs64hi's |x| < |y|/2 holds
		uint64_t q1=(uint64_t)divs64hi((int64_t)(1ULL<<61),(int64_t)ax);
		uint64_t rem1=0-q1*ax;                                         // 2^125 - q1*ax < ax: the low word is all of it
		uint64_t q2=powi_div63(rem1,ax,&sticky);                       // the remainder is the sticky part
		bh=(q1<<1)|(q2>>62);                                           // R in (2^125,2^126), times 4
		bl=q2<<2;
		be=-126-ex;                                                    // 1/|x| = R*2^(-124-ex)
	}
	// left to right over the bits of |n|: square, and multiply by the base where the bit is set
	uint64_t h=bh,l=bl;
	int64_t e=be;
	int32_t bit=62-PD_CLZ((uint64_t)un);                               // the bit below the top one (|n| >= 2 here)
	if(n>0) {
		e=2*e+powi_norm(multu64hi(bh,bh),bh*bh,0,0,&h,&l,&sticky);     // the first square of a 64-bit base is one product
		if((un>>bit)&1) {
			e+=be+powi_mul64(&h,&l,bh,&sticky);
		}
		bit--;
	}
	for(;bit>=0;bit--) {
		e=2*e+powi_sqr(&h,&l,&sticky);
		if((un>>bit)&1) {
			e+=be+(n>0?powi_mul64(&h,&l,bh,&sticky):powi_mul(&h,&l,bh,bl,&sticky));
		}
	}
	// round (h:l)*2^e once: (h>>1)*2^(e+65), anything below that only matters as sticky
	sticky|=(h&1)|l;
	int64_t em=e+129+PSEUDO_DOUBLE_EXP_BIAS;
	if(em>(1LL<<30)) {
		em=1LL<<30;                                                    // far out of range either way: the helper
	} else if(em<-(1LL<<30)) {                                         // reports the overflow or underflow
		em=-(1LL<<30);
	}
	pseudo_double_i r=pseudo_round_helper((signed_pd_internal)(h>>1),(int32_t)em,sticky==0);
	return negative?pdi_neg(r):r;
}

pseudo_double_i pdi_floor(pseudo_double_i x) {
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	if(e<2) {
		return (((signed_pd_internal)x)<0)?0xC000000000000000ULL+PSEUDO_DOUBLE_EXP_BIAS+2:0; // |x| < 1: -1 or 0
	}
	if(e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS) {
		return x;
	}
	uint64_t m=(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-e))-1;
	pseudo_double_i r=x&~m;
	if(r==0x8000000000000000ULL) {
		return 0xC000000000000000ULL+exponent+1;                       // a mantissa of -2^63 is -2^62 one exponent up
	}
	return r+exponent;
}

pseudo_double_i pdi_ceil(pseudo_double_i x) {
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal mantissa=x&EXP_MASK_INV;
	if(e<2) {
		return (((signed_pd_internal)x)>0)?(1ULL<<62)+PSEUDO_DOUBLE_EXP_BIAS+2:0; // |x| < 1: 1 or 0
	}
	if(e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS) {
		return x;
	}
	uint64_t m=(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-e-1))-1;
	signed_pd_internal vr=((mantissa>>1)+m)&~m;
	return pseudo_round_helper(vr,exponent+1,true);
}

pseudo_double_i pdi_round(pseudo_double_i x) {
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal mantissa=x&EXP_MASK_INV;
	if(e<1) {
		return 0;
	}
	if(e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS) {
		return x;
	}
	uint64_t add=1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-e-2);
	uint64_t m=(add<<1)-1;
	signed_pd_internal vr=((mantissa>>1)+(mantissa>0?add:add-1))&~m;
	return pseudo_round_helper(vr,exponent+1,true);
}

pseudo_double_i pdi_trunc(pseudo_double_i x) {
	return (((signed_pd_internal)x)<0)?pdi_ceil(x):pdi_floor(x);
}

// x - trunc(x): the fractional part, with the sign of x (as Rust's f64::fract), so fract(x) + trunc(x) = x exactly
pseudo_double_i pdi_fract(pseudo_double_i x) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;            // x = vx*2^(e-64): the low 64-e bits of vx are the fraction
	if(e<=1) {
		return x;
	}
	if(e>=64-PSEUDO_DOUBLE_EXP_BITS) {
		return 0;                                                      // no mantissa bits below the binary point: an integer
	}
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	signed_pd_internal one=(signed_pd_internal)1<<(64-e);              // 1.0 on vx's scale (64-e <= 62 here)
	signed_pd_internal fr=vx&(one-1);                                  // x - floor(x), in [0,1)
	if(vx<0 && fr!=0) {
		fr-=one;                                                       // x - ceil(x) for negative x, in (-1,0)
	}
	return pseudo_round_helper(fr,(int32_t)(x&EXP_MASK),true);         // exact: just normalises
}

// ---- exact scaling by powers of 10 (used by both decimal input and output) ----

static const uint64_t pd_pow10[20]={1ULL,10ULL,100ULL,1000ULL,10000ULL,100000ULL,1000000ULL,10000000ULL,
	100000000ULL,1000000000ULL,10000000000ULL,100000000000ULL,1000000000000ULL,10000000000000ULL,100000000000000ULL,
	1000000000000000ULL,10000000000000000ULL,100000000000000000ULL,1000000000000000000ULL,10000000000000000000ULL};

// (hi:lo)*2^*f times c (2 <= c <= 10^19): the exact 192-bit product, cut to its top 128 bits (hi normalised to bit 63)
static inline void pd_mul10k(uint64_t* hi, uint64_t* lo, int32_t* f, uint64_t c, uint64_t* sticky) {
	uint64_t p2=multu64hi(*hi,c), p1=*hi*c, p0=*lo*c, t=multu64hi(*lo,c);
	p1+=t;
	p2+=p1<t;
	int32_t z=PD_CLZ(p2);                                              // p2 >= 1, as hi >= 2^63 and c >= 2
	*hi=(p2<<z)|((p1>>1)>>(63-z));                                     // (x>>1)>>(63-z) is x>>(64-z), and 0 for z=0
	*lo=(p1<<z)|((p0>>1)>>(63-z));
	*sticky|=p0<<z;
	*f+=64-z;
}

// (hi:lo)*2^*f divided by c (2 <= c <= 10^9): long division in 32-bit limbs, with two extra limbs for 64 more
// quotient bits, then the top 128 bits of the quotient. The remainder and anything cut off go into *sticky
static inline void pd_div10k(uint64_t* hi, uint64_t* lo, int32_t* f, uint64_t c, uint64_t* sticky) {
	uint64_t a[6]={*hi>>32,*hi&0xFFFFFFFFULL,*lo>>32,*lo&0xFFFFFFFFULL,0,0};
	uint64_t q[6],r=0;
	for(int i=0;i<6;i++) {
		uint64_t cur=(r<<32)|a[i];                                     // r < c < 2^30, so this fits
		q[i]=cur/c;
		r=cur%c;
	}
	uint64_t w2=(q[0]<<32)|q[1], w1=(q[2]<<32)|q[3], w0=(q[4]<<32)|q[5]; // floor((hi:lo)*2^64/c)
	int32_t z=PD_CLZ(w2);                                              // w2 is in [2^33,2^63), so 1 <= z <= 30
	*hi=(w2<<z)|(w1>>(64-z));
	*lo=(w1<<z)|(w0>>(64-z));
	*sticky|=(w0<<z)|r;
	*f-=z;
}

// (hi:lo)*2^*f times 10^k, in steps of 10^19 up and 10^9 down
static inline void pd_scale10(uint64_t* hi, uint64_t* lo, int32_t* f, int32_t k, uint64_t* sticky) {
	for(;k>=19;k-=19) {
		pd_mul10k(hi,lo,f,10000000000000000000ULL,sticky);
	}
	if(k>0) {
		pd_mul10k(hi,lo,f,pd_pow10[k],sticky);
	}
	for(;k<=-9;k+=9) {
		pd_div10k(hi,lo,f,1000000000ULL,sticky);
	}
	if(k<0) {
		pd_div10k(hi,lo,f,pd_pow10[-k],sticky);
	}
}

// value = (mhi:mlo)*10^e, negated if negative, rounded once. sticky means non-zero digits after the mantissa were
// dropped (so the mantissa is slightly below the true value). Rounds to nearest, ties to even when the value is
// exact, like strtod.
static inline pseudo_double_i decimal_to_pdi(uint64_t mhi, uint64_t mlo, bool negative, int32_t e, bool sticky) {
	if(mhi==0 && mlo==0) {
		return 0;
	}
	// value = (hi:lo)*2^f, hi normalised to bit 63
	uint64_t hi,lo;
	int32_t f;
	if(mhi==0) {
		int c=PD_CLZ(mlo);
		hi=mlo<<c;
		lo=0;
		f=-c-64;
	} else {
		int c=PD_CLZ(mhi);
		hi=c==0?mhi:(mhi<<c)|(mlo>>(64-c));
		lo=mlo<<c;
		f=-c;
	}
	// far outside the range either way: clamp, so the scaling below stays short (the helper then reports it)
	const int32_t elim=(int32_t)(((1LL<<(PSEUDO_DOUBLE_EXP_BITS-1))*1233)>>12)+60;
	e=e>elim?elim:(e<-elim?-elim:e);
	uint64_t dropped=sticky;  // non-zero if bits have been dropped, so (hi:lo) is a floor
	pd_scale10(&hi,&lo,&f,e,&dropped);
	dropped|=(hi&1)|lo;       // dropped by the shift below
	// hi>>1 is a positive normalised mantissa; if bits were dropped it is a floor, so half up is correct
	pseudo_double_i r=pseudo_round_helper((signed_pd_internal)(hi>>1),f+129+PSEUDO_DOUBLE_EXP_BIAS,dropped==0);
	return negative?pdi_neg(r):r;
}

pseudo_double_i int64fixed10_to_pdi(int64_t d, int32_t e) {
	return decimal_to_pdi(0,d<0?0-(uint64_t)d:(uint64_t)d,d<0,e,false);
}

// parses [+-]digits[.digits][(e|E)[+-]digits], with at least one digit before any exponent.
// The value is formed as mant*10^exp10 from the first 38 significant digits (in 128 bits) and rounded once, ties to
// even when exact (as strtod does); later digits only set a sticky flag. So the result is correctly rounded except
// within about 10^-38 of a rounding boundary, which needs at least 39 digits to reach.
pseudo_double_i string_to_pdi(const char* str) {
	bool neg=false;
	if(*str=='-' || *str=='+') {
		neg=(*str=='-');
		str++;
	}
	uint64_t mhi=0,mlo=0;                // first 38 significant digits: below 10^38 < 2^127
	int32_t exp10=0;                     // value = mant*10^exp10
	int32_t digits=0;                    // significant digits kept (leading zeros don't count)
	bool any=false;                      // at least one digit seen
	bool sticky=false;                   // a non-zero digit was dropped
	bool point=false;                    // past the decimal point
	for(;;str++) {
		if(*str=='.' && !point) {
			point=true;
			continue;
		}
		if(*str<'0' || *str>'9') {
			break;
		}
		any=true;
		if(digits<38) {
			uint64_t d=(uint64_t)(*str-'0');
			uint64_t l=mlo*10, c=multu64hi(mlo,10);  // mant = mant*10 + d
			l+=d;
			c+=(l<d);
			mhi=mhi*10+c;
			mlo=l;
			digits+=(mhi!=0 || mlo!=0);
			exp10-=point;
		} else {
			exp10+=!point;               // a dropped integer digit still scales the value
			sticky|=(*str!='0');
		}
	}
	if(!any) {
		PD_DO_ERROR_SYNTAX;
		return PD_NAN;
	}
	if(*str=='e' || *str=='E') {
		str++;
		bool eneg=false;
		if(*str=='-' || *str=='+') {
			eneg=(*str=='-');
			str++;
		}
		if(*str<'0' || *str>'9') {
			PD_DO_ERROR_SYNTAX;
			return PD_NAN;
		}
		int32_t e=0;
		for(;*str>='0' && *str<='9';str++) {
			if(e<100000000) {
				e=e*10+(*str-'0');       // saturate: far beyond any exponent range, and no int overflow
			}
		}
		exp10+=eneg?-e:e;
	}
	if(*str!='\0') {
		PD_DO_ERROR_SYNTAX;
		return PD_NAN;
	}
	return decimal_to_pdi(mhi,mlo,neg,exp10,sticky);
}

// ---- conversion to decimal strings ----
// |x| is scaled by a power of 10 to a 128-bit value (a floor) with a sticky flag for every bit dropped, giving 19
// significant digits plus 64 more fraction bits. Rounding to up to 19 digits is then exact (ties to even when the
// value is exact), except within about 2^-59 of a unit in the 19th digit of a rounding boundary. Digits asked for
// beyond the 19th are written as 0

typedef struct {
	uint64_t d19;                                                      // 19 digits, in [10^18,10^19)
	uint64_t frac;                                                     // the next 64 bits, as 0.64
	uint64_t sticky;                                                   // nonzero: more bits below frac
	int32_t exp10;                                                     // power of 10 of the first digit
} pd_decimal_t;

// |x| for a positive nonzero x: the floor of |x|*10^(18-exp10), as d19 then frac and sticky
static inline void pd_to_decimal(pseudo_double_i a, pd_decimal_t* r) {
	int32_t ex=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;           // |x| = ax*2^(ex-64), in [2^(ex-2),2^(ex-1))
	uint64_t ax=(uint64_t)(a&EXP_MASK_INV);
	int32_t d0=(int32_t)(((int64_t)(ex-2)*1292913986LL)>>32);          // floor((ex-2)*log10(2)): exp10 is d0 or d0+1
	for(;;) {
		// |x|*10^(17-d0) is then in [10^17,10^19), below 2^64 (if d0 were off, try again)
		uint64_t hi=ax<<1, lo=0, st=0;
		int32_t f=ex-129;                                              // |x| = (hi:lo)*2^f
		pd_scale10(&hi,&lo,&f,17-d0,&st);
		int32_t s=-f;                                                  // value = (hi:lo)/2^s
		if(s<64) {
			d0++;
			continue;
		}
		uint64_t dd,frac;
		if(s==64) {
			dd=hi;
			frac=lo;
		} else if(s<128) {
			dd=hi>>(s-64);
			frac=(hi<<(128-s))|(lo>>(s-64));
			st|=lo<<(128-s);
		} else {
			d0--;
			continue;
		}
		if(dd>=pd_pow10[19]) {
			d0++;
			continue;
		}
		if(dd<pd_pow10[17]) {
			d0--;
			continue;
		}
		if(dd<pd_pow10[18]) {                                          // 18 digits: one more from frac (exact)
			dd=dd*10+multu64hi(frac,10);
			frac*=10;
			r->exp10=d0;
		} else {
			r->exp10=d0+1;
		}
		r->d19=dd;
		r->frac=frac;
		r->sticky=st;
		return;
	}
}

// d19 rounded to n significant digits (1 <= n <= 19), to nearest with ties to even: returns q in [10^(n-1),10^n),
// with the power of 10 of its first digit in *exp10
static inline uint64_t pd_decimal_round(const pd_decimal_t* r, int n, int32_t* exp10) {
	uint64_t q,up;
	if(n>=19) {
		q=r->d19;
		up=r->frac>(1ULL<<63) || (r->frac==(1ULL<<63) && (r->sticky!=0 || (q&1)));
	} else {
		uint64_t p=pd_pow10[19-n], tail=r->d19%p, half=p>>1;
		q=r->d19/p;
		up=tail>half || (tail==half && ((r->frac|r->sticky)!=0 || (q&1)));
	}
	q+=up;
	*exp10=r->exp10;
	if(q==pd_pow10[n]) {                                               // rounded up through a power of 10
		q=pd_pow10[n-1];
		*exp10+=1;
	}
	return q;
}

// the fewest significant digits that read back (with string_to_pdi) as the positive x. For n digits the candidates
// are the two n-digit decimals either side of x: the one x rounds to if it reads back, else the other one (which
// happens when x is close to the middle, and at powers of 2 where the gap below is half the gap above). Some n
// digits read back for all n from some point on, so a binary search finds it. Precision p bits always reads back
// with 2+floor(p*log10(2)) digits
static inline int pd_decimal_n_reads_back(pseudo_double_i a, const pd_decimal_t* r, int n, uint64_t* q, int32_t* exp10) {
	uint64_t c=pd_decimal_round(r,n,exp10);
	if(decimal_to_pdi(0,c,false,*exp10-n+1,false)==a) {
		*q=c;
		return 1;
	}
	uint64_t f=r->d19/pd_pow10[19-n];                                  // the n digits below x
	uint64_t other=(c==f && *exp10==r->exp10)?f+1:f;                   // whichever of f and f+1 c isn't
	if(decimal_to_pdi(0,other,false,r->exp10-n+1,false)==a) {
		*exp10=r->exp10;
		if(other==pd_pow10[n]) {
			other=pd_pow10[n-1];
			*exp10+=1;
		}
		*q=other;
		return 1;
	}
	return 0;
}

static inline int pd_decimal_shortest(pseudo_double_i a, const pd_decimal_t* r, uint64_t* q, int32_t* exp10) {
	int lo=1, hi=2+(int)(((63-PSEUDO_DOUBLE_EXP_BITS)*1233)>>12);
	while(lo<hi) {
		int mid=(lo+hi)>>1;
		uint64_t c;
		int32_t e;
		if(pd_decimal_n_reads_back(a,r,mid,&c,&e)) {
			hi=mid;
		} else {
			lo=mid+1;
		}
	}
	if(!pd_decimal_n_reads_back(a,r,lo,q,exp10)) {
		*q=pd_decimal_round(r,lo,exp10);                               // (not reached)
	}
	return lo;
}

// the exact decimal digits of a positive integer x below 2^127, into digs (at least 40 chars): returns how many
static int pd_int_digits(pseudo_double_i a, char* digs) {
	int32_t e=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS-64;         // x = m*2^e, e <= 64
	uint64_t m=(uint64_t)(a&EXP_MASK_INV), hi, lo;
	if(e>=64) {
		hi=m<<(e-64);
		lo=0;
	} else if(e>0) {
		hi=m>>(64-e);
		lo=m<<e;
	} else {
		hi=0;
		lo=m>>-e;                                                      // exact: x is an integer
	}
	char t[45];
	int k=0;
	do {
		// divide (hi:lo) by 10^9 in 32-bit limbs, taking the remainder's 9 digits
		uint64_t w[4]={hi>>32,hi&0xFFFFFFFFULL,lo>>32,lo&0xFFFFFFFFULL}, r=0;
		for(int i=0;i<4;i++) {
			uint64_t cur=(r<<32)|w[i];
			w[i]=cur/1000000000ULL;
			r=cur%1000000000ULL;
		}
		hi=(w[0]<<32)|w[1];
		lo=(w[2]<<32)|w[3];
		for(int j=0;j<9;j++) {
			t[k++]=(char)('0'+r%10);
			r/=10;
		}
	} while((hi|lo)!=0);
	while(k>1 && t[k-1]=='0') {
		k--;                                                           // leading zeros
	}
	for(int i=0;i<k;i++) {
		digs[i]=t[k-1-i];
	}
	return k;
}

// whether a positive x is below 2^127 (and so at least 2^precision, an integer once it has 20 digits)
static inline int pd_below_2_127(pseudo_double_i a) {
	return (int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS-64<=64;
}

// ---- exact digits for long outputs (more than 19 significant digits) ----
// |x|*10^k is formed exactly as a big integer (|x| = m*2^e, so m*10^k shifted by e), then written out in decimal.
// This covers anything up to PD_BIG_LIMBS*32 bits; longer still (only for huge exponent widths or precisions)
// falls back to 19 correct digits followed by zeros

#ifndef PD_BIG_LIMBS
#define PD_BIG_LIMBS 160                                               // 5120 bits
#endif

typedef struct {
	uint32_t w[PD_BIG_LIMBS];                                          // little endian
	int n;                                                             // limbs in use
} pd_big_t;

static inline void pd_big_mul(pd_big_t* b, uint32_t m) {
	uint64_t c=0;
	for(int i=0;i<b->n;i++) {
		c+=(uint64_t)b->w[i]*m;
		b->w[i]=(uint32_t)c;
		c>>=32;
	}
	if(c!=0) {
		b->w[b->n++]=(uint32_t)c;
	}
}

static inline uint32_t pd_big_div(pd_big_t* b, uint32_t d) {             // returns the remainder
	uint64_t r=0;
	for(int i=b->n-1;i>=0;i--) {
		uint64_t cur=(r<<32)|b->w[i];
		b->w[i]=(uint32_t)(cur/d);
		r=cur%d;
	}
	while(b->n>0 && b->w[b->n-1]==0) {
		b->n--;
	}
	return (uint32_t)r;
}

// |x|*10^k for a positive x, as decimal digits (most significant first, no leading zeros, "0" for 0) of its floor,
// returning the count (0 if it won't fit). *round_info: 0 exact, 1 below a half, 2 exactly a half, 3 above a half
static int pd_exact_scaled(pseudo_double_i a, int32_t k, char* digs, int maxdigs, int* round_info) {
	int32_t e=(int32_t)(a&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS-64;         // |x| = m*2^e
	uint64_t m=(uint64_t)(a&EXP_MASK_INV);
	if(k<0 || (int64_t)64+(int64_t)k*10/3+1+(e>0?e:0)>(int64_t)PD_BIG_LIMBS*32-32) {
		return 0;
	}
	pd_big_t b;
	b.w[0]=(uint32_t)m;
	b.w[1]=(uint32_t)(m>>32);
	b.n=2;
	for(;k>=9;k-=9) {
		pd_big_mul(&b,1000000000u);
	}
	if(k>0) {
		pd_big_mul(&b,(uint32_t)pd_pow10[k]);
	}
	*round_info=0;
	if(e>0) {                                                          // shift left
		int32_t ws=e>>5, bs=e&31;
		if(bs!=0) {
			pd_big_mul(&b,1u<<bs);
		}
		if(ws>0) {
			for(int i=b.n-1;i>=0;i--) {
				b.w[i+ws]=b.w[i];
			}
			for(int i=0;i<ws;i++) {
				b.w[i]=0;
			}
			b.n+=ws;
		}
	} else if(e<0) {                                                   // shift right, classifying what drops off
		int32_t s=-e;
		int half=0, rest=0;
		for(int32_t i=0;i<s;i++) {
			int32_t limb=i>>5;
			int bit=limb<b.n?(int)((b.w[limb]>>(i&31))&1):0;
			if(i==s-1) {
				half=bit;
			} else {
				rest|=bit;
			}
		}
		*round_info=half?(rest?3:2):(rest?1:0);
		int32_t ws=s>>5, bs=s&31;
		int nn=b.n-ws;
		if(nn<=0) {
			b.n=0;
		} else {
			for(int i=0;i<nn;i++) {
				uint64_t v=b.w[i+ws];
				if(bs!=0 && i+ws+1<b.n) {
					v|=(uint64_t)b.w[i+ws+1]<<32;
				}
				b.w[i]=(uint32_t)(v>>bs);
			}
			b.n=nn;
		}
		while(b.n>0 && b.w[b.n-1]==0) {
			b.n--;
		}
	}
	// to decimal, 9 digits at a time from the bottom
	char t[PD_BIG_LIMBS*10+10];
	int nt=0;
	do {
		uint32_t r=pd_big_div(&b,1000000000u);
		for(int j=0;j<9;j++) {
			t[nt++]=(char)('0'+r%10);
			r/=10;
		}
	} while(b.n>0);
	while(nt>1 && t[nt-1]=='0') {
		nt--;
	}
	if(nt>maxdigs) {
		return 0;
	}
	for(int i=0;i<nt;i++) {
		digs[i]=t[nt-1-i];
	}
	return nt;
}

// rounds the digit string digs[0..len-1] up by one in the last place: returns 1 if it carried out of the top
// (digs is then 1000...)
static inline int pd_digits_inc(char* digs, int len) {
	for(int i=len-1;i>=0;i--) {
		if(digs[i]!='9') {
			digs[i]++;
			return 0;
		}
		digs[i]='0';
	}
	digs[0]='1';
	return 1;
}

// |x| to exactly n significant digits (n > 19), correctly rounded: digs gets n digits, *ex the power of 10 of the
// first one (exp10 is the first digit's power before rounding). Returns 0 if it won't fit
static int pd_exact_sig(pseudo_double_i a, int n, int32_t exp10, char* digs, int maxdigs, int32_t* ex) {
	int32_t k=n-1-exp10;
	int info,len;
	if(k>=0) {
		len=pd_exact_scaled(a,k,digs,maxdigs,&info);                   // exactly n digits
		if(len!=n) {
			return 0;
		}
	} else {
		// at least 10^20, so an integer: all its digits, then round at n from the ones after
		len=pd_exact_scaled(a,0,digs,maxdigs,&info);
		if(len!=exp10+1) {
			return 0;
		}
		int rest=0;
		for(int i=n+1;i<len;i++) {
			rest|=digs[i]!='0';
		}
		info=digs[n]>'5'?3:(digs[n]=='5'?(rest?3:2):((digs[n]!='0' || rest)?1:0));
	}
	*ex=exp10;
	if(info==3 || (info==2 && ((digs[n-1]-'0')&1))) {
		*ex+=pd_digits_inc(digs,n);
	}
	return 1;
}

// ---- formatting ----

typedef struct {
	char* buf;
	size_t size;
	size_t len;                                                        // what the whole result needs
} pd_out_t;

static inline void pd_out(pd_out_t* o, char c) {
	if(o->len+1<o->size) {
		o->buf[o->len]=c;
	}
	o->len++;
}

// the significant digits digs[0..n-1] (first at 10^x), as %f with p digits after the point; digits past n are 0
static void pd_out_fixed(pd_out_t* o, const char* digs, int n, int32_t x, int32_t p, int alt, char point) {
	if(x<0) {
		pd_out(o,'0');
	} else {
		for(int32_t i=0;i<=x;i++) {
			pd_out(o,i<n?digs[i]:'0');
		}
	}
	if(p>0 || alt) {
		pd_out(o,point);
	}
	for(int32_t j=1;j<=p;j++) {
		int32_t i=x+j;
		pd_out(o,(i>=0 && i<n)?digs[i]:'0');
	}
}

// the same as %e with p digits after the point
static void pd_out_exp(pd_out_t* o, const char* digs, int n, int32_t x, int32_t p, int alt, int upper, char point) {
	pd_out(o,digs[0]);
	if(p>0 || alt) {
		pd_out(o,point);
	}
	for(int32_t i=1;i<=p;i++) {
		pd_out(o,i<n?digs[i]:'0');
	}
	pd_out(o,upper?'E':'e');
	pd_out(o,x<0?'-':'+');
	uint32_t ux=x<0?0u-(uint32_t)x:(uint32_t)x;
	char t[12];
	int k=0;
	do {
		t[k++]=(char)('0'+ux%10);
		ux/=10;
	} while(ux!=0);
	if(k<2) {
		t[k++]='0';                                                    // at least two exponent digits, as printf
	}
	while(k>0) {
		pd_out(o,t[--k]);
	}
}

// the digits of q, exactly n of them (n <= 19)
static inline void pd_digits(char* digs, uint64_t q, int n) {
	for(int i=n-1;i>=0;i--) {
		digs[i]=(char)('0'+q%10);
		q/=10;
	}
}

// Formats x as printf's %.{precision}{conv} does, for conv 'e', 'f', 'g' (or 'E', 'F', 'G'), with PD_FMT_ flags,
// into buf (size bytes, NUL terminated when size > 0). A negative precision means the shortest digits that read
// back as x: 'e' and 'f' write all of them, and 'g' takes whichever of the two is shorter (fixed on a tie), as
// C++'s std::to_chars does. point is the decimal point character. Returns the length the whole result needs, as
// snprintf does. Correctly rounded (ties to even when exact) to any number of digits, except that outputs needing
// more than PD_BIG_LIMBS*32 bits of working (only possible with very wide exponents or huge precisions) get 19
// correct digits and then zeros
int pdi_format_point(char* buf, size_t size, pseudo_double_i x, char conv, int precision, int flags, char point) {
	pd_out_t o={buf,size,0};
	int upper=(conv>='A' && conv<='Z');
	char c=(char)(conv|0x20);
	int alt=(flags&PD_FMT_ALT)!=0;
	int negative=x!=PD_NAN && ((signed_pd_internal)x)<0;
	if(negative) {
		pd_out(&o,'-');
	} else if(flags&PD_FMT_PLUS) {
		pd_out(&o,'+');
	} else if(flags&PD_FMT_SPACE) {
		pd_out(&o,' ');
	}
	if(x==PD_NAN) {
		const char* s=upper?"NAN":"nan";
		while(*s) {
			pd_out(&o,*s++);
		}
	} else {
		pseudo_double_i a=negative?pdi_neg(x):x;
		char digs[PD_BIG_LIMBS*10+10];
		int n=1;
		int32_t ex=0;                                                  // power of 10 of digs[0]
		pd_decimal_t r={0,0,0,0};
		if(x==0) {
			digs[0]='0';
		} else {
			pd_to_decimal(a,&r);
		}
		if(precision<0) {
			// shortest digits that read back as x
			if(x!=0) {
				uint64_t q;
				n=pd_decimal_shortest(a,&r,&q,&ex);
				pd_digits(digs,q,n);
			}
			int32_t fp=n-1-ex>0?n-1-ex:0;                              // digits needed after the point
			if(c=='g') {
				int32_t flen=(ex>=0?ex+1:1)+(fp>0?fp+1:0);
				int32_t elen=n+(n>1)+4;                                // digits, point, e, sign, two exponent digits
				for(int32_t ae=(ex<0?-ex:ex)/100;ae!=0;ae/=10) {
					elen++;                                            // more exponent digits
				}
				c=elen<flen?'e':'f';
			}
			if(c=='e') {
				pd_out_exp(&o,digs,n,ex,n-1,alt,upper,point);
			} else {
				if(ex+1>n && pd_below_2_127(a)) {
					n=pd_int_digits(a,digs);                           // an integer: its exact digits, as std::to_chars
					ex=n-1;
				}
				pd_out_fixed(&o,digs,n,ex,fp,alt,point);
			}
		} else if(c=='e') {
			int want=precision+1;
			if(x!=0) {
				n=want<19?want:19;
				if(want<=19 || !pd_exact_sig(a,n=want,r.exp10,digs,(int)sizeof digs,&ex)) {
					n=want<19?want:19;
					pd_digits(digs,pd_decimal_round(&r,n,&ex),n);
				}
			}
			pd_out_exp(&o,digs,n,ex,precision,alt,upper,point);
		} else if(c=='f') {
			if(x!=0) {
				int64_t want=(int64_t)r.exp10+1+precision;             // significant digits down to 10^-precision
				if(want>19 && want<(int64_t)sizeof digs && pd_exact_sig(a,(int)want,r.exp10,digs,(int)sizeof digs,&ex)) {
					n=(int)want;                                       // more than 19 digits: exact
				} else if(want>19 && r.exp10>=19 && pd_below_2_127(a)) {
					n=pd_int_digits(a,digs);                           // an integer of 20 to 39 digits: exact
					ex=n-1;
				} else if(want>=1) {
					n=want<19?(int)want:19;
					pd_digits(digs,pd_decimal_round(&r,n,&ex),n);
				} else {
					// below the last place: rounds to 0 or to one unit of 10^-precision (ties to even, so to 0)
					int up=want==0 && (r.d19>5*pd_pow10[18] || (r.d19==5*pd_pow10[18] && (r.frac|r.sticky)!=0));
					digs[0]=up?'1':'0';
					ex=up?-precision:0;
				}
			}
			pd_out_fixed(&o,digs,n,ex,precision,alt,point);
		} else {
			// 'g': precision significant digits (at least 1), then %e if the exponent is below -4 or at least the
			// precision, else %f; trailing zeros dropped unless '#'
			int pr=precision==0?1:precision;
			if(x!=0) {
				if(pr<=19 || !pd_exact_sig(a,n=pr,r.exp10,digs,(int)sizeof digs,&ex)) {
					n=pr<19?pr:19;
					pd_digits(digs,pd_decimal_round(&r,n,&ex),n);
				}
			}
			int last=n-1;                                              // last nonzero digit
			while(last>0 && digs[last]=='0') {
				last--;
			}
			if(ex<-4 || ex>=pr) {
				pd_out_exp(&o,digs,n,ex,alt?pr-1:last,alt,upper,point);
			} else {
				int32_t fp=alt?pr-1-ex:(last-ex>0?last-ex:0);
				pd_out_fixed(&o,digs,n,ex,fp,alt,point);
			}
		}
	}
	if(size>0) {
		buf[o.len<size?o.len:size-1]='\0';
	}
	return (int)o.len;
}

int pdi_format(char* buf, size_t size, pseudo_double_i x, char conv, int precision, int flags) {
	return pdi_format_point(buf,size,x,conv,precision,flags,'.');
}

// the shortest string that reads back as x (string_to_pdi), fixed or scientific, whichever is shorter. 32 bytes
// is always enough
int pdi_to_string(char* buf, size_t size, pseudo_double_i x) {
	return pdi_format_point(buf,size,x,'g',-1,0,'.');
}

#if !defined(PD_HW_DIV128) && PSEUDO_DOUBLE_EXP_BITS>=8
// pdi_div's rare path: the remainder u*2^64-a*d is in [0,8d), so step a up to the exact floor
pseudo_double_i pdi_div_exact(unsigned_pd_internal a, unsigned_pd_internal u, unsigned_pd_internal d, unsigned_pd_internal sg, int32_t exp_max) {
	unsigned_pd_internal rl=0-a*d;
	unsigned_pd_internal rh=u-multu64hi(a,d)-(rl!=0);
	while(rh!=0 || rl>=d) {
		rh-=(rl<d);
		rl-=d;
		a++;
	}
	unsigned_pd_internal q=(a^sg)-sg-(sg&(rl!=0));
	return pseudo_round_helper((signed_pd_internal)q,exp_max,false);
}
#endif
