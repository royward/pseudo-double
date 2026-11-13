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

#ifndef PSEUDO_DOUBLE_H
#define PSEUDO_DOUBLE_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>

#ifdef _MSC_VER
#include <intrin.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

// Can set any of these or just use the defaults
//#define PSEUDO_DOUBLE_EXP_BITS 10
// This will default to return PD_NAN
//#define PD_DO_ERROR_OVERFLOW return PD_NAN
//#define PD_DO_ERROR_UNDERFLOW
//#define PD_DO_ERROR_RANGE return PD_NAN

// use this to turn off all the intrinsics and use the code in pseudo_double_portable.h
//#define PD_FORCE_PORTABLE 1

// Mantissa (two's complement, normalised on its magnitude, so x and -x have the same exponent)
// [.1000..01 .. .1100] = (-0.5 .. -0.25]
// [.0100 .. .0111] = [0.25 .. 0.5)
// .1000..00 (-0.5) is never a valid mantissa, so negation is exact and never changes the exponent

#define PSEUDO_DOUBLE_TOTAL_BITS 64

#ifndef PSEUDO_DOUBLE_EXP_BITS
// Recommend 8,10,11 or 16 for this. Other values should work but are untested. Less than 8 will mean transcendental functions are inaccurate in the last few bits.
#define PSEUDO_DOUBLE_EXP_BITS 10
#endif

#define EXP_MASK ((1LL<<PSEUDO_DOUBLE_EXP_BITS)-1)
#define EXP_MASK_INV (~((1ULL<<PSEUDO_DOUBLE_EXP_BITS)-1))

typedef uint64_t pseudo_double_i;
typedef uint64_t unsigned_pd_internal;
typedef int64_t signed_pd_internal;

#if defined(PD_FORCE_PORTABLE)
// Any C99 compiler
#include "pseudo_double_portable.h"
#define PD_CLZ(x) pd_clz64_portable(x)
static inline unsigned_pd_internal multu64hi(unsigned_pd_internal x,unsigned_pd_internal y) {return pd_umulh_portable(x,y);}
static inline signed_pd_internal mults64hi(signed_pd_internal x,signed_pd_internal y) {return pd_mulh_portable(x,y);}
#elif defined(_MSC_VER) && defined(_M_X64)
static inline int PD_CLZ(unsigned __int64 x) {unsigned long i; _BitScanReverse64(&i,x); return 63-(int)i;}
#define PD_HW_DIV128 1
static inline unsigned_pd_internal multu64hi(unsigned_pd_internal x,unsigned_pd_internal y) {unsigned __int64 hi; _umul128(x,y,&hi); return hi;}
static inline signed_pd_internal mults64hi(signed_pd_internal x,signed_pd_internal y) {__int64 hi; _mul128(x,y,&hi); return hi;}
static inline signed_pd_internal divs64hi(signed_pd_internal x,signed_pd_internal y) {__int64 r; __int64 q=_div128(x,0,y,&r);return q-((r!=0)&((r^y)<0));}
#elif defined(_MSC_VER) && defined(_M_ARM64)
#include "pseudo_double_portable.h"
#define PD_CLZ(x) ((int)_CountLeadingZeros64(x))
static inline unsigned_pd_internal multu64hi(unsigned_pd_internal x,unsigned_pd_internal y) {return __umulh(x,y);}
static inline signed_pd_internal mults64hi(signed_pd_internal x,signed_pd_internal y) {return __mulh(x,y);}
#elif defined(__SIZEOF_INT128__)
#define PD_CLZ(x) __builtin_clzll(x)
__extension__ typedef __int128 signed_large_pd_internal;
__extension__ typedef unsigned __int128 unsigned_large_pd_internal;
static inline signed_pd_internal mults64hi(signed_pd_internal x,signed_pd_internal y) {return (signed_pd_internal)(((signed_large_pd_internal)x*y)>>64);}
static inline unsigned_pd_internal multu64hi(unsigned_pd_internal x,unsigned_pd_internal y) {return (unsigned_pd_internal)(((unsigned_large_pd_internal)x*y)>>64);}
#if defined(__x86_64__)
#define PD_HW_DIV128 1
// floor(x*2^64/y) with the hardware divide; remainder comes free in rdx.
// Requires abs(x)<abs(y)/2 or (abs(x)==abs(y)/2,x and y different signs) (guaranteed by the shift right by 2 in pdi_div) so the quotient fits, otherwise idiv faults.
static inline signed_pd_internal divs64hi(signed_pd_internal x,signed_pd_internal y) {
	signed_pd_internal q,r;
	__asm__("idivq %4" : "=a"(q), "=d"(r) : "a"((signed_pd_internal)0), "d"(x), "rm"(y) : "cc");
	return q-((r!=0)&((r^y)<0));
}
#endif
#else
#include "pseudo_double_portable.h"
#define PD_CLZ(x) pd_clz64_portable(x)
static inline unsigned_pd_internal multu64hi(unsigned_pd_internal x,unsigned_pd_internal y) {return pd_umulh_portable(x,y);}
static inline signed_pd_internal mults64hi(signed_pd_internal x,signed_pd_internal y) {return pd_mulh_portable(x,y);}
#endif
// high 64 bits of signed x times unsigned y (the unsigned product, corrected for x's sign); same for every platform
static inline signed_pd_internal mulsu64hi(signed_pd_internal x,unsigned_pd_internal y) {return (signed_pd_internal)(multu64hi((unsigned_pd_internal)x,y)-(x<0?y:0));}

#ifndef PD_HW_DIV128
// floor(x*2^64/y) for targets with no 128/64 divide (ARM64, RISC-V, ...), where x/y on 128 bits calls a
// general 128/128 routine (__divti3). Same contract as the idiv version: |x|<|y|/2, or |x|==|y|/2 with opposite signs.
// Uses no division: a reciprocal from a table plus Newton steps, then a 2-by-1 step simplified for the zero low word
// (Moller & Granlund, "Improved division by invariant integers", IEEE Trans. Computers 60(2), 2011, Algorithms 2 and 4).
static const uint16_t pd_recip_table[256]={ // floor((2^19-3*2^8)/d9) for d9=256..511
	2045,2037,2029,2021,2013,2005,1998,1990,1983,1975,1968,1960,1953,1946,1938,1931,
	1924,1917,1910,1903,1896,1889,1883,1876,1869,1863,1856,1849,1843,1836,1830,1824,
	1817,1811,1805,1799,1792,1786,1780,1774,1768,1762,1756,1750,1745,1739,1733,1727,
	1722,1716,1710,1705,1699,1694,1688,1683,1677,1672,1667,1661,1656,1651,1646,1641,
	1636,1630,1625,1620,1615,1610,1605,1600,1596,1591,1586,1581,1576,1572,1567,1562,
	1558,1553,1548,1544,1539,1535,1530,1526,1521,1517,1513,1508,1504,1500,1495,1491,
	1487,1483,1478,1474,1470,1466,1462,1458,1454,1450,1446,1442,1438,1434,1430,1426,
	1422,1418,1414,1411,1407,1403,1399,1396,1392,1388,1384,1381,1377,1374,1370,1366,
	1363,1359,1356,1352,1349,1345,1342,1338,1335,1332,1328,1325,1322,1318,1315,1312,
	1308,1305,1302,1299,1295,1292,1289,1286,1283,1280,1276,1273,1270,1267,1264,1261,
	1258,1255,1252,1249,1246,1243,1240,1237,1234,1231,1228,1226,1223,1220,1217,1214,
	1211,1209,1206,1203,1200,1197,1195,1192,1189,1187,1184,1181,1179,1176,1173,1171,
	1168,1165,1163,1160,1158,1155,1153,1150,1148,1145,1143,1140,1138,1135,1133,1130,
	1128,1125,1123,1121,1118,1116,1113,1111,1109,1106,1104,1102,1099,1097,1095,1092,
	1090,1088,1086,1083,1081,1079,1077,1074,1072,1070,1068,1066,1064,1061,1059,1057,
	1055,1053,1051,1049,1047,1044,1042,1040,1038,1036,1034,1032,1030,1028,1026,1024,};
static inline signed_pd_internal divs64hi(signed_pd_internal x,signed_pd_internal y) {
	unsigned_pd_internal ux=x<0?0-(unsigned_pd_internal)x:(unsigned_pd_internal)x;
	unsigned_pd_internal uy=y<0?0-(unsigned_pd_internal)y:(unsigned_pd_internal)y;
	int s=PD_CLZ(uy);
	unsigned_pd_internal d=uy<<s, u=ux<<s;                                // d normalised, u<=d/2
	// v=floor((2^128-1)/d)-2^64
	unsigned_pd_internal d0=d&1, d40=(d>>24)+1, d63=(d>>1)+d0;
	unsigned_pd_internal v0=pd_recip_table[(d>>55)&255];
	unsigned_pd_internal v1=(v0<<11)-((v0*v0*d40)>>40)-1;
	unsigned_pd_internal v2=(v1<<13)+((v1*((1ULL<<60)-v1*d40))>>47);
	unsigned_pd_internal e=((v2>>1)&(0-d0))-v2*d63;
	unsigned_pd_internal v3=(multu64hi(v2,e)>>1)+(v2<<31);
	unsigned_pd_internal lo=v3*d+d;
	unsigned_pd_internal v=v3-(multu64hi(v3,d)+(lo<d)+d);
	// q=floor(u*2^64/d): the estimate is exact or one too big (about 90% of the time, so the fix is branch-free)
	unsigned_pd_internal q0=v*u, q=multu64hi(v,u)+u+1;
	unsigned_pd_internal r=0-q*d;                                         // u*2^64-q*d, mod 2^64
	unsigned_pd_internal m=0-(unsigned_pd_internal)(r>q0);
	q+=m;
	r+=m&d;
	if(r>=d) {q++; r-=d;}                                                 // practically never
	unsigned_pd_internal sg=0-(((unsigned_pd_internal)(x^y))>>63);        // branch-free sign: signs are often random
	return (signed_pd_internal)((q^sg)-sg-(sg&(r!=0)));
}
#endif
#if defined(__GNUC__) || defined(__clang__)
#define PD_UNLIKELY(x) __builtin_expect(!!(x),0)
#else
#define PD_UNLIKELY(x) (x)
#endif

#define PSEUDO_DOUBLE_EXP_BIAS (1<<(PSEUDO_DOUBLE_EXP_BITS-1))
#define PSEUDO_DOUBLE_HALF_ULP ((1ULL<<(PSEUDO_DOUBLE_EXP_BITS-1))-1)
#define PD_NAN ((pseudo_double_i)(-1))

#ifndef PD_ERROR_CHECK
// Setting this to 0 will turn off most overflow/underflow/range checking and result in a tiny speed increase
// Probably not worth it in most scenarios
#define PD_ERROR_CHECK 1
#endif
#ifndef PD_DO_ERROR_OVERFLOW
#define PD_DO_ERROR_OVERFLOW return PD_NAN
#endif
#ifndef PD_DO_ERROR_UNDERFLOW
#define PD_DO_ERROR_UNDERFLOW
#endif
#ifndef PD_DO_ERROR_INTEGER_OVERFLOW
#define PD_DO_ERROR_INTEGER_OVERFLOW
#endif
#ifndef PD_DO_ERROR_RANGE
#define PD_DO_ERROR_RANGE return PD_NAN
#endif
#ifndef PD_DO_ERROR_SYNTAX
#define PD_DO_ERROR_SYNTAX return PD_NAN
#endif

static inline signed_pd_internal shift_left_signed(signed_pd_internal x, int shift) {
	if(shift>=0) {
		return (signed_pd_internal)((unsigned_pd_internal)x<<shift);
	}
	return x>>-shift;
}

static inline unsigned_pd_internal shift_left_unsigned(unsigned_pd_internal x, int shift) {
	if(shift>=0) {
		return x<<shift;
	}
	return x>>-shift;
}

static inline pseudo_double_i pdi_neg(pseudo_double_i x) {
	return (pseudo_double_i)(-(signed_pd_internal)(x&EXP_MASK_INV))+(x&EXP_MASK);
}

static inline pseudo_double_i pdi_abs(pseudo_double_i x) {
	signed_pd_internal m=(signed_pd_internal)(x&EXP_MASK_INV);
	signed_pd_internal s=m>>(PSEUDO_DOUBLE_TOTAL_BITS-1);
	return (pseudo_double_i)((m^s)-s)+(x&EXP_MASK);
}

static inline pseudo_double_i pdi_max(pseudo_double_i x, pseudo_double_i y) {
	int neg=((unsigned_pd_internal)y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1);
	if((x^y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return neg?x:y;
	}
	// signs are the same, check exponent
	int expdiff=(signed_pd_internal)((x&EXP_MASK)-(y&EXP_MASK));
	if(expdiff!=0) {
		return ((expdiff>0)^neg)?x:y;
	}
	// exponents are the same so don't need to mask off, check mantissa
	return (x>y)?x:y;
}

static inline pseudo_double_i pdi_min(pseudo_double_i x, pseudo_double_i y) {
	int neg=((unsigned_pd_internal)y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1);
	if((x^y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return neg?y:x;
	}
	// signs are the same, check exponent
	int expdiff=(signed_pd_internal)((x&EXP_MASK)-(y&EXP_MASK));
	if(expdiff!=0) {
		return ((expdiff>0)^neg)?y:x;
	}
	// exponents are the same so don't need to mask off, check mantissa
	return (x>y)?y:x;
}

static inline bool pdi_gt(pseudo_double_i x, pseudo_double_i y) {
	bool neg=((unsigned_pd_internal)y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1);
	if((x^y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return neg;
	}
	// signs are the same, check exponent
	int expdiff=(signed_pd_internal)((x&EXP_MASK)-(y&EXP_MASK));
	if(expdiff!=0) {
		return (expdiff>0)^neg;
	}
	// exponents are the same so don't need to mask off, check mantissa
	return (x>y);
}

static inline bool pdi_gte(pseudo_double_i x, pseudo_double_i y) {
	bool neg=((unsigned_pd_internal)y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1);
	if((x^y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return neg;
	}
	// signs are the same, check exponent
	int expdiff=(signed_pd_internal)((x&EXP_MASK)-(y&EXP_MASK));
	if(expdiff!=0) {
		return (expdiff>0)^neg;
	}
	// exponents are the same so don't need to mask off, check mantissa
	return (x>=y);
}

static inline pseudo_double_i pseudo_round_helper(signed_pd_internal vr,int32_t exp_max, bool ties_even) {
	if(vr==0) {
		// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
		return (pseudo_double_i)0;
	}
	int32_t leading_bits=PD_CLZ(vr>0?vr:~vr)-1;
	int32_t new_exponent=exp_max-leading_bits;
	vr=(signed_pd_internal)((unsigned_pd_internal)vr<<leading_bits);
	signed_pd_internal vr2=(((unsigned_pd_internal)vr)+(PSEUDO_DOUBLE_HALF_ULP+(ties_even?((vr>>PSEUDO_DOUBLE_EXP_BITS)&1):1)))&EXP_MASK_INV;
	if(PD_UNLIKELY(vr2==(signed_pd_internal)(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-1)))) {
		// +2^63 (rounded up past the top) or -2^63 (a negative power of two, which the leading bit count shifts one
		// too far, or rounded down to one): halve it, keeping vr's sign
		vr2=(signed_pd_internal)((((unsigned_pd_internal)vr2)>>1)|(((unsigned_pd_internal)vr)&(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-1))));
		new_exponent+=1;
	}
	if(PD_UNLIKELY((uint32_t)new_exponent>(uint32_t)EXP_MASK)) { // one test for both ends of the range
		if(new_exponent<0) {
#if PD_ERROR_CHECK
			PD_DO_ERROR_UNDERFLOW;
#endif
			return (pseudo_double_i)0;
		}
#if PD_ERROR_CHECK
		PD_DO_ERROR_OVERFLOW;
#endif
	}
	return (pseudo_double_i)((vr2&EXP_MASK_INV)+new_exponent);
}

static inline pseudo_double_i pdi_add(pseudo_double_i x, pseudo_double_i y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	int32_t ydiffx=expy-expx;
	if(ydiffx>=(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return y;
	}
	if(ydiffx<=-(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return x;
	}
	signed_pd_internal vx=((signed_pd_internal)(x&EXP_MASK_INV))>>1;
	signed_pd_internal vy=((signed_pd_internal)(y&EXP_MASK_INV))>>1;
	bool y_big=ydiffx>=0;                                              // selects rather than branches: GCC turns these into cmovs
	signed_pd_internal big=y_big?vy:vx, small=y_big?vx:vy;
	int32_t shift=y_big?ydiffx:-ydiffx;                                // 0..62
	int32_t exp_max=(y_big?expy:expx)+1;
	unsigned_pd_internal lost=((unsigned_pd_internal)small<<1)<<(63-shift); // the bits shifted out (0 for a shift of 0)
	return pseudo_round_helper(big+(small>>shift),exp_max,lost==0);   // bits lost: floor, so half up; exact: ties to even
}

static inline pseudo_double_i pdi_sub(pseudo_double_i x, pseudo_double_i y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	int32_t ydiffx=expy-expx;
	if(ydiffx>=(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return pdi_neg(y);
	}
	if(ydiffx<=-(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return x;
	}
	signed_pd_internal vx=((signed_pd_internal)(x&EXP_MASK_INV))>>1;
	signed_pd_internal vy=-(((signed_pd_internal)(y&EXP_MASK_INV))>>1);
	bool y_big=ydiffx>=0;
	signed_pd_internal big=y_big?vy:vx, small=y_big?vx:vy;
	int32_t shift=y_big?ydiffx:-ydiffx;
	int32_t exp_max=(y_big?expy:expx)+1;
	unsigned_pd_internal lost=((unsigned_pd_internal)small<<1)<<(63-shift);
	return pseudo_round_helper(big+(small>>shift),exp_max,lost==0);
}

static inline pseudo_double_i pdi_mult(pseudo_double_i x, pseudo_double_i y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	signed_pd_internal vy=(signed_pd_internal)(y&EXP_MASK_INV);
	signed_pd_internal vr=mults64hi(vx,vy);
	return pseudo_round_helper(vr,expx+expy-PSEUDO_DOUBLE_EXP_BIAS,(uint64_t)vx*(uint64_t)vy==0);
}

#if !defined(PD_HW_DIV128) && PSEUDO_DOUBLE_EXP_BITS>=8
// Division without a 128/64 divide. A 64-bit divide gives a reciprocal good to about 31 bits, and one Goldschmidt
// step takes the quotient to within 8 below the exact one. That rounds the same way as the exact quotient unless a
// rounding boundary lies in between, which the low bits show; only then (1-5% of calls) is the exact quotient
// worked out. Needs at least 8 exponent bits.
// Rare path of pdi_div (in pseudo_double.c): a<=u*2^64/d<a+8, made exact, signed by sg and rounded
pseudo_double_i pdi_div_exact(unsigned_pd_internal a, unsigned_pd_internal u, unsigned_pd_internal d, unsigned_pd_internal sg, int32_t exp_max);
static inline pseudo_double_i pdi_div(pseudo_double_i x, pseudo_double_i y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	signed_pd_internal vy=(signed_pd_internal)(y&EXP_MASK_INV);
	if(vy==0) { // leave this one in to avoid division by zero signal
		PD_DO_ERROR_RANGE;
		return PD_NAN;
	}
	int32_t exp_max=expx-expy+PSEUDO_DOUBLE_EXP_BIAS+2;
	unsigned_pd_internal ux=vx<0?0-(unsigned_pd_internal)vx:(unsigned_pd_internal)vx;   // |vx|<2^63 (|vx|/4 is the dividend)
	unsigned_pd_internal uy=vy<0?0-(unsigned_pd_internal)vy:(unsigned_pd_internal)vy;   // in [2^62,2^63)
	unsigned_pd_internal d=uy<<1;                       // D=d/2^63 in [1,2)
	unsigned_pd_internal u2=ux;                         // 2u, where the exact quotient is T=u*2^64/d<2^63
	// R=floor((2^64-1)/((d>>32)+1))/2^33 <= 1/D; eps=1-D*R < 1.5*2^-31, as 0.64 rounded down
	unsigned_pd_internal r=((~(unsigned_pd_internal)0)/((d>>32)+1))<<31;
	unsigned_pd_internal ph=multu64hi(d,r), pl=d*r;
	unsigned_pd_internal e=(((1ULL<<63)-ph-(pl!=0))<<1)|((0-pl)>>63);
	unsigned_pd_internal a=multu64hi(u2,r);             // T*(1-eps), rounded down
	a+=multu64hi(a,e);                                  // T*(1-eps^2), rounded down: T-7.01 < a <= T
	unsigned_pd_internal sg=0-(((unsigned_pd_internal)(vx^vy))>>63);
	unsigned_pd_internal q=(a^sg)-sg;                   // exact floor is in [q,q+7] (q>=0) or [q-8,q] (q<0)
	// The rounding step for this quotient: pseudo_round_helper drops EXP_BITS-lb bits, where lb is 0, 1 or 2 as
	// T>=2^62, T>=2^61 or less. Found from the operands, so off the critical path; near the power of two where lb
	// changes, both sides round to that power of two, so the choice there does not matter.
	int lb=(int)(u2<(d>>1))+(int)(u2<(d>>2));
	unsigned_pd_internal cell=(1ULL<<PSEUDO_DOUBLE_EXP_BITS)>>lb;
	// distance from a to the next rounding boundary, on the side the exact quotient lies (above a, so below -a)
	unsigned_pd_internal w=((cell>>1)-1-a-sg)&(cell-1);
	if(w<8) {
		return pdi_div_exact(a,u2>>1,d,sg,exp_max);
	}
	return pseudo_round_helper((signed_pd_internal)q,exp_max,false);
}
#else
static inline pseudo_double_i pdi_div(pseudo_double_i x, pseudo_double_i y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	signed_pd_internal vy=(signed_pd_internal)(y&EXP_MASK_INV);
	if(vy==0) { // leave this one in to avoid division by zero signal
		PD_DO_ERROR_RANGE;
		return PD_NAN;
	}
	signed_pd_internal vr=divs64hi(vx>>2,vy);
	return pseudo_round_helper(vr,expx-expy+PSEUDO_DOUBLE_EXP_BIAS+2,false);
}
#endif

static inline pseudo_double_i pdi_ldexp(pseudo_double_i x, int y) {
	if(x==0) {
		return 0;
	}
	int32_t expx=x&EXP_MASK;
#if PD_ERROR_CHECK
	if(expx+y>(int32_t)EXP_MASK) {
		PD_DO_ERROR_OVERFLOW;
	}
#endif
	if(expx+y<0) {
#if PD_ERROR_CHECK
		PD_DO_ERROR_UNDERFLOW;
#endif
		return (pseudo_double_i)0;
	}
	return x+y;
}

pseudo_double_i double_to_pdi(double d);
pseudo_double_i int64_to_pdi(int64_t d);
pseudo_double_i uint64_to_pdi(uint64_t d);
double pdi_to_double(pseudo_double_i d);
int64_t pdi_to_int64(pseudo_double_i d);
uint64_t pdi_to_uint64(pseudo_double_i d);

pseudo_double_i int64fixed10_to_pdi(int64_t d, int32_t e);
pseudo_double_i int64fixed2_to_pdi(int64_t d, int32_t e);
int64_t pdi_to_int64fixed2(pseudo_double_i d, int32_t e);

pseudo_double_i pdi_floor(pseudo_double_i x);
pseudo_double_i pdi_ceil(pseudo_double_i x);
pseudo_double_i pdi_round(pseudo_double_i x);
pseudo_double_i pdi_trunc(pseudo_double_i x);
pseudo_double_i pdi_sqrt(pseudo_double_i x);
pseudo_double_i pdi_inv_sqrt(pseudo_double_i x);
pseudo_double_i pdi_exp2(pseudo_double_i x);
pseudo_double_i pdi_exp10(pseudo_double_i x);
pseudo_double_i pdi_exp(pseudo_double_i x);
pseudo_double_i pdi_log2(pseudo_double_i x);
pseudo_double_i pdi_log(pseudo_double_i x);
pseudo_double_i pdi_log10(pseudo_double_i x);
pseudo_double_i pdi_pow(pseudo_double_i x, pseudo_double_i y);
pseudo_double_i pdi_sin_rev(pseudo_double_i x);
pseudo_double_i pdi_cos_rev(pseudo_double_i x);
pseudo_double_i pdi_tan_rev(pseudo_double_i x);
pseudo_double_i pdi_atan2_rev(pseudo_double_i y, pseudo_double_i x);
pseudo_double_i pdi_sin(pseudo_double_i x);
pseudo_double_i pdi_cos(pseudo_double_i x);
pseudo_double_i pdi_tan(pseudo_double_i x);
pseudo_double_i pdi_atan2(pseudo_double_i y, pseudo_double_i x);
pseudo_double_i pdi_expm1(pseudo_double_i x);
pseudo_double_i pdi_ln_1p(pseudo_double_i x);
pseudo_double_i pdi_powi(pseudo_double_i x, int32_t n);
void pdi_sin_cos_rev(pseudo_double_i x, pseudo_double_i* sin_out, pseudo_double_i* cos_out);
void pdi_sin_cos(pseudo_double_i x, pseudo_double_i* sin_out, pseudo_double_i* cos_out);
pseudo_double_i pdi_atan_rev(pseudo_double_i x);
pseudo_double_i pdi_atan(pseudo_double_i x);
pseudo_double_i pdi_asin_rev(pseudo_double_i x);
pseudo_double_i pdi_acos_rev(pseudo_double_i x);
pseudo_double_i pdi_asin(pseudo_double_i x);
pseudo_double_i pdi_acos(pseudo_double_i x);
pseudo_double_i pdi_sinh(pseudo_double_i x);
pseudo_double_i pdi_cosh(pseudo_double_i x);
pseudo_double_i pdi_tanh(pseudo_double_i x);
pseudo_double_i pdi_asinh(pseudo_double_i x);
pseudo_double_i pdi_acosh(pseudo_double_i x);
pseudo_double_i pdi_atanh(pseudo_double_i x);
pseudo_double_i pdi_cbrt(pseudo_double_i x);
pseudo_double_i pdi_fract(pseudo_double_i x);
pseudo_double_i pdi_hypot(pseudo_double_i x, pseudo_double_i y);
pseudo_double_i pdi_hypot3(pseudo_double_i x, pseudo_double_i y, pseudo_double_i z);
pseudo_double_i pdi_inv_hypot(pseudo_double_i x, pseudo_double_i y);
pseudo_double_i pdi_inv_hypot3(pseudo_double_i x, pseudo_double_i y, pseudo_double_i z);
pseudo_double_i pdi_mul_add(pseudo_double_i a, pseudo_double_i b, pseudo_double_i c); // a*b+c
pseudo_double_i pdi_mul_sub(pseudo_double_i a, pseudo_double_i b, pseudo_double_i c); // a*b-c
pseudo_double_i pdi_sub_mul(pseudo_double_i a, pseudo_double_i b, pseudo_double_i c); // c-a*b
pseudo_double_i string_to_pdi(const char* str);

// x is a 2.62 unsigned fixed in the range (1,4)
// result is 1.63 unsigned fixed in the range (0.5,1)
uint64_t inv_sqrt64_fixed(uint64_t x);

// x is a 0.64 unsigned fixed in the range [0,1)
// result is 2.62 unsigned fixed in the range [1,2)
uint64_t exp2_64_fixed(uint64_t x);

// x is a 0.64 unsigned fixed in the range [0,1)
// calculate log2(x+1)
// result is 0.64 unsigned fixed in the range [0,1)
uint64_t log2_64_fixed(uint64_t x);

// x is a 2.62 unsigned fixed in the range [0,1]
// calculate sin with x in quarter turns
// result is 2.62 unsigned fixed in the range [0,1]
uint64_t sin_rev_64_fixed(uint64_t x);

// x is a 2.62 unsigned fixed in the range [0,1]
// calculate atan with x in eighth turns
// result is 2.62 unsigned fixed in the range [0,1]
uint64_t atan_rev_64_fixed(uint64_t x);

// x is a 2.62 unsigned fixed in the range [0,1] (quarter turns)
// sets sin and cos of x quarter turns, as 2.62 unsigned fixed in the range [0,1]
void sin_cos_rev_64_fixed(uint64_t x, uint64_t* sin_out, uint64_t* cos_out);

// x is a 3.61 unsigned fixed in the range [1,8)
// result is x^(-1/3) as 1.63 unsigned fixed in the range (0.5,1]
uint64_t inv_cbrt64_fixed(uint64_t x);

#define PD_FMT_PLUS 1    // '+': always write a sign
#define PD_FMT_SPACE 2   // ' ': a space where a plus sign would go
#define PD_FMT_ALT 4     // '#': always write the decimal point, and %g keeps its trailing zeros

// Formats x as printf's %.{precision}{conv} does, for conv 'e', 'f', 'g' (or 'E', 'F', 'G'), with PD_FMT_ flags,
// into buf (size bytes, NUL terminated when size > 0). A negative precision means the shortest digits that read back
// as x ('g' then takes fixed or scientific, whichever is shorter, as std::to_chars does). Returns the length the
// whole result needs, as snprintf does. Correctly rounded, ties to even when exact
int pdi_format(char* buf, size_t size, pseudo_double_i x, char conv, int precision, int flags);
// the same, with point as the decimal point character
int pdi_format_point(char* buf, size_t size, pseudo_double_i x, char conv, int precision, int flags, char point);
// the shortest string that reads back as x with string_to_pdi (32 bytes is always enough)
int pdi_to_string(char* buf, size_t size, pseudo_double_i x);

/* =========================================================================================================================================================================
 * below this line is the C wrapper
 * =========================================================================================================================================================================
 */

typedef struct {
	pseudo_double_i val;
} pseudo_double;

static inline pseudo_double create_pseudo_double_from_internal(pseudo_double_i x) {
	pseudo_double ret;
	ret.val=x;
	return ret;
}

static inline pseudo_double pd_neg(pseudo_double x) {return create_pseudo_double_from_internal(pdi_neg(x.val));}
static inline pseudo_double pd_abs(pseudo_double x) {return create_pseudo_double_from_internal(pdi_abs(x.val));}
static inline bool pd_gt(pseudo_double x, pseudo_double y) {return pdi_gt(x.val,y.val);}
static inline bool pd_gte(pseudo_double x, pseudo_double y) {return pdi_gte(x.val,y.val);}
static inline pseudo_double pd_max(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_max(x.val,y.val));}
static inline pseudo_double pd_min(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_min(x.val,y.val));}
static inline pseudo_double pd_sub(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_sub(x.val,y.val));}
static inline pseudo_double pd_add(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_add(x.val,y.val));}
static inline pseudo_double pd_mult(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_mult(x.val,y.val));}
static inline pseudo_double pd_div(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_div(x.val,y.val));}
static inline pseudo_double pd_ldexp(pseudo_double x, int y) {return create_pseudo_double_from_internal(pdi_ldexp(x.val,y));}
static inline pseudo_double double_to_pd(double d) {return create_pseudo_double_from_internal(double_to_pdi(d));}
static inline pseudo_double int64_to_pd(int64_t d) {return create_pseudo_double_from_internal(int64_to_pdi(d));}
static inline pseudo_double uint64_to_pd(uint64_t d) {return create_pseudo_double_from_internal(uint64_to_pdi(d));}
static inline double pd_to_double(pseudo_double d) {return pdi_to_double(d.val);}
static inline int64_t pd_to_int64(pseudo_double d) {return pdi_to_int64(d.val);}
static inline uint64_t pd_to_uint64(pseudo_double d) {return pdi_to_uint64(d.val);}
static inline pseudo_double int64fixed10_to_pd(int64_t d, int32_t e) {return create_pseudo_double_from_internal(int64fixed10_to_pdi(d,e));}
static inline pseudo_double int64fixed2_to_pd(int64_t d, int32_t e) {return create_pseudo_double_from_internal(int64fixed2_to_pdi(d,e));}
static inline int64_t pd_to_int64fixed2(pseudo_double d, int32_t e) {return pdi_to_int64fixed2(d.val,e);}
static inline pseudo_double pd_floor(pseudo_double x) {return create_pseudo_double_from_internal(pdi_floor(x.val));}
static inline pseudo_double pd_ceil(pseudo_double x) {return create_pseudo_double_from_internal(pdi_ceil(x.val));}
static inline pseudo_double pd_round(pseudo_double x) {return create_pseudo_double_from_internal(pdi_round(x.val));}
static inline pseudo_double pd_trunc(pseudo_double x) {return create_pseudo_double_from_internal(pdi_trunc(x.val));}
static inline pseudo_double pd_sqrt(pseudo_double x) {return create_pseudo_double_from_internal(pdi_sqrt(x.val));}
static inline pseudo_double pd_inv_sqrt(pseudo_double x) {return create_pseudo_double_from_internal(pdi_inv_sqrt(x.val));}
static inline pseudo_double pd_exp2(pseudo_double x) {return create_pseudo_double_from_internal(pdi_exp2(x.val));}
static inline pseudo_double pd_exp10(pseudo_double x) {return create_pseudo_double_from_internal(pdi_exp10(x.val));}
static inline pseudo_double pd_exp(pseudo_double x) {return create_pseudo_double_from_internal(pdi_exp(x.val));}
static inline pseudo_double pd_log2(pseudo_double x) {return create_pseudo_double_from_internal(pdi_log2(x.val));}
static inline pseudo_double pd_log(pseudo_double x) {return create_pseudo_double_from_internal(pdi_log(x.val));}
static inline pseudo_double pd_log10(pseudo_double x) {return create_pseudo_double_from_internal(pdi_log10(x.val));}
static inline pseudo_double pd_pow(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_pow(x.val,y.val));}
static inline pseudo_double pd_sin_rev(pseudo_double x) {return create_pseudo_double_from_internal(pdi_sin_rev(x.val));}
static inline pseudo_double pd_cos_rev(pseudo_double x) {return create_pseudo_double_from_internal(pdi_cos_rev(x.val));}
static inline pseudo_double pd_tan_rev(pseudo_double x) {return create_pseudo_double_from_internal(pdi_tan_rev(x.val));}
static inline pseudo_double pd_atan2_rev(pseudo_double y, pseudo_double x) {return create_pseudo_double_from_internal(pdi_atan2_rev(y.val,x.val));}
static inline pseudo_double pd_sin(pseudo_double x) {return create_pseudo_double_from_internal(pdi_sin(x.val));}
static inline pseudo_double pd_cos(pseudo_double x) {return create_pseudo_double_from_internal(pdi_cos(x.val));}
static inline pseudo_double pd_tan(pseudo_double x) {return create_pseudo_double_from_internal(pdi_tan(x.val));}
static inline pseudo_double pd_atan2(pseudo_double y, pseudo_double x) {return create_pseudo_double_from_internal(pdi_atan2(y.val,x.val));}
static inline bool pd_gt_zero(pseudo_double x) {return (((signed_pd_internal)x.val)>0);}
static inline bool pd_gte_zero(pseudo_double x) {return (((signed_pd_internal)x.val)>=0);}
static inline bool pd_lt_zero(pseudo_double x) {return (((signed_pd_internal)x.val)<0);}
static inline bool pd_lte_zero(pseudo_double x) {return (((signed_pd_internal)x.val)<=0);}
static inline bool pd_eq_zero(pseudo_double x) {return (x.val==0);}
static inline bool pd_neq_zero(pseudo_double x) {return (x.val!=0);}
static inline pseudo_double string_to_pd(const char* str) {return create_pseudo_double_from_internal(string_to_pdi(str));}
static inline pseudo_double pd_expm1(pseudo_double x) {return create_pseudo_double_from_internal(pdi_expm1(x.val));}
static inline pseudo_double pd_ln_1p(pseudo_double x) {return create_pseudo_double_from_internal(pdi_ln_1p(x.val));}
static inline pseudo_double pd_powi(pseudo_double x, int32_t n) {return create_pseudo_double_from_internal(pdi_powi(x.val,n));}
static inline void pd_sin_cos_rev(pseudo_double x, pseudo_double* sin_out, pseudo_double* cos_out) {pdi_sin_cos_rev(x.val,&sin_out->val,&cos_out->val);}
static inline void pd_sin_cos(pseudo_double x, pseudo_double* sin_out, pseudo_double* cos_out) {pdi_sin_cos(x.val,&sin_out->val,&cos_out->val);}
static inline pseudo_double pd_atan_rev(pseudo_double x) {return create_pseudo_double_from_internal(pdi_atan_rev(x.val));}
static inline pseudo_double pd_atan(pseudo_double x) {return create_pseudo_double_from_internal(pdi_atan(x.val));}
static inline pseudo_double pd_asin_rev(pseudo_double x) {return create_pseudo_double_from_internal(pdi_asin_rev(x.val));}
static inline pseudo_double pd_acos_rev(pseudo_double x) {return create_pseudo_double_from_internal(pdi_acos_rev(x.val));}
static inline pseudo_double pd_asin(pseudo_double x) {return create_pseudo_double_from_internal(pdi_asin(x.val));}
static inline pseudo_double pd_acos(pseudo_double x) {return create_pseudo_double_from_internal(pdi_acos(x.val));}
static inline pseudo_double pd_sinh(pseudo_double x) {return create_pseudo_double_from_internal(pdi_sinh(x.val));}
static inline pseudo_double pd_cosh(pseudo_double x) {return create_pseudo_double_from_internal(pdi_cosh(x.val));}
static inline pseudo_double pd_tanh(pseudo_double x) {return create_pseudo_double_from_internal(pdi_tanh(x.val));}
static inline pseudo_double pd_asinh(pseudo_double x) {return create_pseudo_double_from_internal(pdi_asinh(x.val));}
static inline pseudo_double pd_acosh(pseudo_double x) {return create_pseudo_double_from_internal(pdi_acosh(x.val));}
static inline pseudo_double pd_atanh(pseudo_double x) {return create_pseudo_double_from_internal(pdi_atanh(x.val));}
static inline pseudo_double pd_cbrt(pseudo_double x) {return create_pseudo_double_from_internal(pdi_cbrt(x.val));}
static inline pseudo_double pd_fract(pseudo_double x) {return create_pseudo_double_from_internal(pdi_fract(x.val));}
static inline pseudo_double pd_hypot(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_hypot(x.val,y.val));}
static inline pseudo_double pd_hypot3(pseudo_double x, pseudo_double y, pseudo_double z) {return create_pseudo_double_from_internal(pdi_hypot3(x.val,y.val,z.val));}
static inline pseudo_double pd_inv_hypot(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_inv_hypot(x.val,y.val));}
static inline pseudo_double pd_inv_hypot3(pseudo_double x, pseudo_double y, pseudo_double z) {return create_pseudo_double_from_internal(pdi_inv_hypot3(x.val,y.val,z.val));}
static inline pseudo_double pd_mul_add(pseudo_double a, pseudo_double b, pseudo_double c) {return create_pseudo_double_from_internal(pdi_mul_add(a.val,b.val,c.val));}
static inline pseudo_double pd_mul_sub(pseudo_double a, pseudo_double b, pseudo_double c) {return create_pseudo_double_from_internal(pdi_mul_sub(a.val,b.val,c.val));}
static inline pseudo_double pd_sub_mul(pseudo_double a, pseudo_double b, pseudo_double c) {return create_pseudo_double_from_internal(pdi_sub_mul(a.val,b.val,c.val));}
static inline int pd_to_string(char* buf, size_t size, pseudo_double x) {return pdi_to_string(buf,size,x.val);}
static inline int pd_format(char* buf, size_t size, pseudo_double x, char conv, int precision, int flags) {return pdi_format(buf,size,x.val,conv,precision,flags);}

#ifdef __cplusplus
}
#endif

#endif // PSEUDO_DOUBLE_H
