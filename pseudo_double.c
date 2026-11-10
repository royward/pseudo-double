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
	uint64_t m=(vx<0)?0-(uint64_t)vx:(uint64_t)vx;      // magnitude in [2^62,2^63]; 0x8000... is just 2^63
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
// END generated by gen_exp2_table.py - do not edit

// x is a 0.64 unsigned fixed in the range [0,1)
// result is 2.62 unsigned fixed in the range [1,2)
uint64_t exp2_64_fixed(uint64_t x) {
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
	return (t+(multu64hi(t,q)>>5))>>1;              // (t + t*q), 1.63 -> 2.62
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

// 2^t, where t = (h:l)*2^(sh-128) is a 128-bit signed value (h:l) with |h:l| in [2^124,2^126]
// shared by pdi_exp, pdi_exp10 and pdi_pow, which all form t in 128 bits so the only rounding is at the end
static inline pseudo_double_i pdi_exp2_wide(int64_t h, uint64_t l, int32_t sh) {
	if(sh<=-61) {
		return (1ULL<<62)+PSEUDO_DOUBLE_EXP_BIAS+2;                     // |t| < 2^-63: 2^t rounds to 1
	}
	if(sh>=PSEUDO_DOUBLE_EXP_BITS+4) {                                  // |t| >= 2^EXP_BITS: certain overflow or underflow
		if(h<0) {
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
	// split t = ipart + fraction, with ipart = floor(t) and fraction in [0,1) as 0.64
	int32_t ipart;
	uint64_t fraction;
	if(sh>=1) {
		ipart=(int32_t)(h>>(64-sh));
		fraction=((uint64_t)h<<sh)|(l>>(64-sh));
	} else {
		ipart=(int32_t)(h>>63);                                         // 0 or -1
		fraction=(uint64_t)(h>>-sh);
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

// log2(x) for x > 0, before rounding: the result is vr*2^(exp_max-PSEUDO_DOUBLE_EXP_BIAS-64), the form pseudo_round_helper takes
// vr is not necessarily normalised, and is 0 for x = 1
static inline int32_t pdi_log2_unrounded(pseudo_double_i x, signed_pd_internal* vr_out) {
	int32_t e=(int32_t)(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS-2;      // x = (1+f)*2^e
	uint64_t f=(uint64_t)(x&EXP_MASK_INV)<<2;                       // f as 0.64
	signed_pd_internal vr;
	int32_t exp_max;
	if((e==0 && f<(1ULL<<58)) || (e==-1 && f>=0-(1ULL<<58))) {
		// x in [1-2^-7,1+2^-6): the result is small, so it is done separately to keep full relative precision
		int64_t r=(int64_t)(e==0?f<<5:f<<4);                        // r = x-1, scaled by 2^69 (exact)
		if(r==0) {
			*vr_out=0;                                              // log(1) = 0
			return 0;
		}
		// log2(1+r) = r*b0 + r*u*tail, with r normalised (rn) so that the result carries full precision
		int32_t s=PD_CLZ(r>0?r:~r)-1;
		int64_t rn=(int64_t)((uint64_t)r<<s);
		vr=mults64hi(rn,LOG2_C0)+mults64hi(mults64hi(rn,r),log2_tail(r,mults64hi(r,r),5));
		exp_max=PSEUDO_DOUBLE_EXP_BIAS-3-s;
	} else {
		// result = e + log2(1+f): put the integer e above the 0.64 fraction, using as few bits for e as possible
		uint64_t lo=log2_64_fixed_inline(f);
		uint32_t a=(uint32_t)(e<0?~e:e);
		int32_t nb=63-PD_CLZ(((uint64_t)a<<1)|1);                   // bits needed for e, excluding the sign bit
		if(nb==0 && (int64_t)(lo^(uint64_t)(int64_t)e)>=0) {
			vr=(signed_pd_internal)lo;                              // e is 0 or -1 and |result| < 1/2: lo as signed 0.64 is exactly the result
			exp_max=PSEUDO_DOUBLE_EXP_BIAS;
		} else {
			vr=(signed_pd_internal)(((uint64_t)(int64_t)e<<(63-nb))|(lo>>(nb+1)));
			exp_max=PSEUDO_DOUBLE_EXP_BIAS+1+nb;
		}
	}
	*vr_out=vr;
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

#if __GNUG__
// static const pseudo_double_i pdi_inv_log_2_e =int64fixed10_to_pdi(6931471805599453094,-19);
// static const pseudo_double_i pdi_inv_log_2_10=int64fixed10_to_pdi(3010299956639811952,-19);
static const pseudo_double_i pdi_1_div_tau=pdi_div(uint64_to_pdi(1),int64fixed10_to_pdi(6283185307179586477,-18));
static const pseudo_double_i pdi_tau=int64fixed10_to_pdi(6283185307179586477,-18);
#endif


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
	int32_t es,ec;
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
	int32_t es,ec;
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
	// now 0 <= y <= x, and the result is add_const +- atan(r) for r = y/x: a = y and b = x as mantissas in [2^62,2^63)
	int32_t ea=(int32_t)(y&EXP_MASK), eb=(int32_t)(x&EXP_MASK);
	uint64_t ma=y==0?0:(uint64_t)(y&EXP_MASK_INV), mb=(uint64_t)(x&EXP_MASK_INV);
	uint64_t t8=0;                                                     // atan(r) in eighth turns, 1.63
	uint32_t ki=0;
	int64_t q=0;                                                       // for r < 1/64: r = q*2^(ea-eb-62)
	if(ma!=0 && ea-eb>=-7) {
		uint64_t ma2=ma>>(eb-ea);                                      // a aligned to b's exponent (exact: low bits are 0)
		uint32_t a26=(uint32_t)(ma2>>37), b26=(uint32_t)(mb>>37);     // top 26 bits: a26*32 + b26 still fits in 32 bits
		ki=(a26*32+(b26>>1)+2)/b26;                                    // round(32*r); the +2 keeps r >= 1/64 - tiny when ki >= 1
		if(ki!=0) {
			int64_t n=(int64_t)(32*ma2-ki*mb);                         // 32*(a - b*k/32), exact: |n| < mb/2, so the bits above drop out
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
	pseudo_double_i r=pseudo_round_helper(vr,exp_max,false);          // round the magnitude, so atan2(-y,x) = -atan2(y,x)
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

pseudo_double_i pdi_floor(pseudo_double_i x) {
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	if(e<2) {
		return (((signed_pd_internal)x)<0)?0x8000000000000000ULL+PSEUDO_DOUBLE_EXP_BIAS+1:0;
	}
	if(e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS) {
		return x;
	}
	uint64_t m=(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-e))-1;
	return (x&~m)+exponent;
}

pseudo_double_i pdi_ceil(pseudo_double_i x) {
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal mantissa=x&EXP_MASK_INV;
	if(e<2) {
		if(e==1 && (uint64_t)mantissa==0x8000000000000000ULL) { // special test for ceil(-1)=-1
			return x;
		}
		return (((signed_pd_internal)x)>0)?(1ULL<<62)+PSEUDO_DOUBLE_EXP_BIAS+2:0;
	}
	if(e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS) {
		return x;
	}
	uint64_t m=(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-e-1))-1;
	signed_pd_internal vr=((mantissa>>1)+m)&~m;
	int32_t leading_bits=PD_CLZ(vr>0?vr:~vr)-1;
	vr=(signed_pd_internal)((uint64_t)vr<<leading_bits);
	int32_t new_exponent=exponent+1-leading_bits;
#if PD_ERROR_CHECK
	if(new_exponent>EXP_MASK) {
		PD_DO_ERROR_OVERFLOW;
	}
#endif
	return (pseudo_double_i)(vr+new_exponent);
}

pseudo_double_i pdi_round(pseudo_double_i x) {
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal mantissa=x&EXP_MASK_INV;
	if(e<1) {
		if(e==0 && (uint64_t)mantissa==0x8000000000000000ULL) { // special test for round(-0.5)=-1
			return x+1; // hack to multiply by 2
		}
		return 0;
	}
	if(e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS) {
		return x;
	}
	uint64_t add=1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-e-2);
	uint64_t m=(add<<1)-1;
	signed_pd_internal vr=((mantissa>>1)+(mantissa>0?add:add-1))&~m;
	int32_t leading_bits=PD_CLZ(vr>0?vr:~vr)-1;
	vr=(signed_pd_internal)((uint64_t)vr<<leading_bits);
	int32_t new_exponent=exponent+1-leading_bits;
#if PD_ERROR_CHECK
	if(new_exponent>EXP_MASK) {
		PD_DO_ERROR_OVERFLOW;
	}
#endif
	return (pseudo_double_i)(vr+new_exponent);
}

pseudo_double_i pdi_trunc(pseudo_double_i x) {
	return (((signed_pd_internal)x)<0)?pdi_ceil(x):pdi_floor(x);
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
