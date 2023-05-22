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

// Can set any of these or just use the defaults
//#define PSEUDO_DOUBLE_TOTAL_BITS 64
//#define PSEUDO_DOUBLE_EXP_BITS 16
// This will default to 2^(PSEUDO_DOUBLE_EXP_BITS-1)
//#define PSEUDO_DOUBLE_EXP_BIAS 32768
// This will default to return PF_NAN
//#define PF_DO_ERROR_OVERFLOW return PF_NAN
//#define PF_DO_ERROR_UNDERFLOW return 0
//#define PF_DO_ERROR_RANGE return PF_NAN

// Mantissa
// [.1000 .. .1011] = [-0.5 .. -0.25)
// [.0100 .. .0111] = [0.25 .. 0.5)

#define PSEUDO_DOUBLE_TOTAL_BITS 64

#ifndef PSEUDO_DOUBLE_EXP_BITS
// Recommend 8,16 or 32 for this. Other values will work but be less efficient
#define PSEUDO_DOUBLE_EXP_BITS 16
#endif

#define EXP_MASK ((1LL<<PSEUDO_DOUBLE_EXP_BITS)-1)
#define EXP_MASK_INV (~((1ULL<<PSEUDO_DOUBLE_EXP_BITS)-1))
// leaving this conditional in in case we finish the 32 bit case. Currently only works with 64 bits
#if PSEUDO_DOUBLE_TOTAL_BITS==64
typedef uint64_t pseudo_double_i;
typedef uint64_t unsigned_pd_internal;
typedef int64_t signed_pd_internal;
typedef __int128 signed_large_pd_internal;
#define clz __builtin_clzll
#elif PSEUDO_DOUBLE_TOTAL_BITS==32
typedef uint32_t pseudo_double_i;
typedef uint32_t unsigned_pd_internal;
typedef int32_t signed_pd_internal;
typedef int64_t signed_large_pd_internal;
#define clz __builtin_clz
#else
#error PSEUDO_DOUBLE_TOTAL_BITS must be 32 or 64 bits
#endif
#ifndef PSEUDO_DOUBLE_EXP_BIAS
#define PSEUDO_DOUBLE_EXP_BIAS (1U<<(PSEUDO_DOUBLE_EXP_BITS-1))
#endif
#define PSEUDO_DOUBLE_HALF_ULP ((1ULL<<(PSEUDO_DOUBLE_EXP_BITS-1))-1)
#define PF_NAN ((pseudo_double_i)(-1))

#ifndef PF_ERROR_CHECK
// Setting this to 0 will turn off most overflow/underflow/range checking and result in a tiny speed increase
// Probably not worth it in most scenarios
#define PF_ERROR_CHECK 1
#endif
#ifndef PF_DO_ERROR_OVERFLOW
#define PF_DO_ERROR_OVERFLOW return PF_NAN
#endif
#ifndef PF_DO_ERROR_UNDERFLOW
#define PF_DO_ERROR_UNDERFLOW return 0
#endif
#ifndef PF_DO_ERROR_RANGE
#define PF_DO_ERROR_RANGE return PF_NAN
#endif

inline signed_pd_internal shift_left_signed(signed_pd_internal x, int shift) {
	if(shift>=0) {
		return x<<shift;
	}
	return x>>-shift;
}

inline unsigned_pd_internal shift_left_unsigned(unsigned_pd_internal x, int shift) {
	if(shift>=0) {
		return x<<shift;
	}
	return x>>-shift;
}

inline pseudo_double_i pdi_neg(pseudo_double_i x) {
	// check for special cases due to the representation range of two's complement being asymmetric:
	// mantissa of 1000000000... can't be negated directly and will need an increase in the exponent
	// mantissa of 0100000000... has exponent decreased after negation
	// everything else negates as expected
	uint32_t exponent=x&EXP_MASK;
	unsigned_pd_internal mantissa=x&EXP_MASK_INV;
	if((mantissa<<2)==0) {
		// get the high order byte
		uint32_t hi_byte=x>>(PSEUDO_DOUBLE_TOTAL_BITS-8);
		if(hi_byte==0x80) {
#if PF_ERROR_CHECK
			if(exponent==EXP_MASK) {
				PF_DO_ERROR_OVERFLOW;
			}
#endif
			return (mantissa>>1)+exponent+1;
		}
		if(hi_byte==0x40) {
#if PF_ERROR_CHECK
			if(exponent==0) {
				PF_DO_ERROR_UNDERFLOW;
			}
#endif
			return (mantissa<<1)+exponent-1;
		}
	}
	return (-(x&EXP_MASK_INV))+exponent;
}

inline pseudo_double_i pdi_abs(pseudo_double_i x) {
	if(((signed_pd_internal)x)>=0) {
		return x;
	}
	uint32_t exponent=x&EXP_MASK;
	unsigned_pd_internal mantissa=x&EXP_MASK_INV;
	if((mantissa<<2)==0) {
		// get the high order byte
		uint32_t hi_byte=x>>(PSEUDO_DOUBLE_TOTAL_BITS-8);
		if(hi_byte==0x80) {
#if PF_ERROR_CHECK
			if(exponent==EXP_MASK) {
				PF_DO_ERROR_OVERFLOW;
			}
#endif
			return (mantissa>>1)+exponent+1;
		}
	}
	return (-(x&EXP_MASK_INV))+exponent;
}

inline int pdi_gt(pseudo_double_i x, pseudo_double_i y) {
	int neg=((unsigned_pd_internal)y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1);
	if((x^y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return neg;
	}
	// signs are the same, check exponent
	int expdiiff=(signed_pd_internal)((x&EXP_MASK)-(y&EXP_MASK));
	if(expdiiff!=0) {
		return (expdiiff>0)^neg;
	}
	// exponents are the same so don't need to mask off, check mantissa
	return (x>y);
}

inline int pdi_gte(pseudo_double_i x, pseudo_double_i y) {
	int neg=((unsigned_pd_internal)y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1);
	if((x^y)>>(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return neg;
	}
	// signs are the same, check exponent
	int expdiiff=(signed_pd_internal)((x&EXP_MASK)-(y&EXP_MASK));
	if(expdiiff!=0) {
		return (expdiiff>0)^neg;
	}
	// exponents are the same so don't need to mask off, check mantissa
	return (x>=y);
}

inline pseudo_double_i pdi_sub(pseudo_double_i x, pseudo_double_i y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	int32_t ydiffx=expy-expx;
	if(ydiffx>=(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return pdi_neg(y);
	}
	if(ydiffx<=-(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return x;
	}
	int32_t exp_max;
	signed_pd_internal vx=((signed_pd_internal)(x&EXP_MASK_INV))>>1;
	signed_pd_internal vy=((signed_pd_internal)(y&EXP_MASK_INV))>>1;
	if(ydiffx>=0) {
		exp_max=expy;
		vx>>=ydiffx;
	} else {
		exp_max=expx;
		vy>>=-ydiffx;
	}
	exp_max+=1;
	signed_pd_internal vr=(vx-vy+PSEUDO_DOUBLE_HALF_ULP)&~PSEUDO_DOUBLE_HALF_ULP;
	if(vr==0) {
		// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
		return (pseudo_double_i)0;
	}
	int32_t leading_bits=clz(vr>0?vr:~vr)-1;
	if(leading_bits>exp_max) {
		leading_bits=exp_max;
	}
	vr<<=leading_bits;
	int32_t new_exponent=exp_max-leading_bits;
#if PF_ERROR_CHECK
	if(new_exponent>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(new_exponent<0) {
		PF_DO_ERROR_UNDERFLOW;
	}
#endif
	return (pseudo_double_i)((vr&EXP_MASK_INV)+new_exponent);
}

inline pseudo_double_i pdi_add(pseudo_double_i x, pseudo_double_i y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	int32_t ydiffx=expy-expx;
	if(ydiffx>=(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return y;
	}
	if(ydiffx<=-(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
		return x;
	}
	int32_t exp_max;
	signed_pd_internal vx=((signed_pd_internal)(x&EXP_MASK_INV))>>1;
	signed_pd_internal vy=((signed_pd_internal)(y&EXP_MASK_INV))>>1;
	if(ydiffx>=0) {
		exp_max=expy;
		vx>>=ydiffx;
	} else {
		exp_max=expx;
		vy>>=-ydiffx;
	}
	exp_max+=1;
	signed_pd_internal vr=(vx+vy+PSEUDO_DOUBLE_HALF_ULP)&~PSEUDO_DOUBLE_HALF_ULP;
	if(vr==0) {
		// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
		return (pseudo_double_i)0;
	}
	int32_t leading_bits=clz(vr>0?vr:~vr)-1;
	if(leading_bits>exp_max) {
		leading_bits=exp_max;
	}
	vr<<=leading_bits;
	int32_t new_exponent=exp_max-leading_bits;
#if PF_ERROR_CHECK
	if(new_exponent>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(new_exponent<0) {
		PF_DO_ERROR_UNDERFLOW;
	}
#endif
	return (pseudo_double_i)((vr&EXP_MASK_INV)+new_exponent);
}

inline pseudo_double_i pdi_mult(pseudo_double_i x, pseudo_double_i y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	signed_pd_internal vy=(signed_pd_internal)(y&EXP_MASK_INV);
	signed_pd_internal vr=(((signed_large_pd_internal)vx)*vy)>>64;
	if(vr==0) {
		// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
		return (pseudo_double_i)0;
	}
	int32_t leading_bits=clz(vr>0?vr:~vr)-1;
	vr<<=leading_bits;
	int32_t new_exponent=expx+expy-PSEUDO_DOUBLE_EXP_BIAS-leading_bits;
#if PF_ERROR_CHECK
	if(new_exponent>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(new_exponent<0) {
		PF_DO_ERROR_UNDERFLOW;
	}
#endif
	return (pseudo_double_i)((vr&EXP_MASK_INV)+new_exponent);
}

inline pseudo_double_i pdi_div(pseudo_double_i x, pseudo_double_i y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	signed_pd_internal vy=(signed_pd_internal)(y&EXP_MASK_INV);
	if(vy==0) { // leave this one in to avoid division bby zero signal
		PF_DO_ERROR_RANGE;
	}
	signed_pd_internal vr=(((((signed_large_pd_internal)vx)>>2)<<64)/vy);
	if(vr==0) {
		// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
		return (pseudo_double_i)0;
	}
	int32_t leading_bits=clz(vr>0?vr:~vr)-1;
	vr<<=leading_bits;
	int32_t new_exponent=2+expx-expy+PSEUDO_DOUBLE_EXP_BIAS-leading_bits;
#if PF_ERROR_CHECK
	if(new_exponent>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(new_exponent<0) {
		PF_DO_ERROR_UNDERFLOW;
	}
#endif
	return (pseudo_double_i)((vr&EXP_MASK_INV)+new_exponent);
}

inline pseudo_double_i pdi_ldexp(pseudo_double_i x, int y) {
#if PF_ERROR_CHECK
	int32_t expx=x&EXP_MASK;
	if(expx+y>(int32_t)EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(expx+y<0) {
		PF_DO_ERROR_UNDERFLOW;
	}
#endif
	return x+y;
}

pseudo_double_i double_to_pdi(double d);
pseudo_double_i int64_to_pdi(int64_t d);
pseudo_double_i uint64_to_pdi(uint64_t d);
double pdi_to_double(pseudo_double_i d);
int64_t pdi_to_int64(pseudo_double_i d);
uint64_t pdi_to_uint64(pseudo_double_i d);

inline pseudo_double_i int64fixed10_to_pdi(int64_t d, int32_t e) {
	if(d==0) {
		return 0;
	}
	int negative=(d<0);
	int32_t nexp=0;
	while(e>0) {
		int lead_bits=clz(negative?~d:d);
		if(lead_bits<5) {
			// check that there is no overflow
			d>>=(5-lead_bits);
			nexp+=(5-lead_bits);
		}
		d*=10;
		e--;
	}
	while(e<0) {
		int lead_bits=clz(negative?~d:d);
		if(lead_bits>1) {
			// make the number as accurate as possible
			d<<=(lead_bits-1);
			nexp-=(lead_bits-1);
		}
		d/=10;
		e++;
	}
	int lead_bits=clz(d);
	int exp=nexp+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits;
	return ((shift_left_signed(d,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))&EXP_MASK_INV)+exp;
}

pseudo_double_i int64fixed2_to_pdi(int64_t d, int32_t e);
int64_t pdi_to_int64fixed2(pseudo_double_i d, int32_t e);

pseudo_double_i pdi_floor(pseudo_double_i x);
pseudo_double_i pdi_ceil(pseudo_double_i x);
pseudo_double_i pdi_round(pseudo_double_i x);
pseudo_double_i pdi_sqrt(pseudo_double_i x);
pseudo_double_i pdi_inv_sqrt(pseudo_double_i x);
pseudo_double_i pdi_exp2(pseudo_double_i x);
pseudo_double_i pdi_exp(pseudo_double_i x);
pseudo_double_i pdi_log2(pseudo_double_i x);
pseudo_double_i pdi_log(pseudo_double_i x);
pseudo_double_i pdi_log10(pseudo_double_i x);
pseudo_double_i pdi_pow(pseudo_double_i x, pseudo_double_i y);
pseudo_double_i pdi_sin_rev(pseudo_double_i x);
pseudo_double_i pdi_cos_rev(pseudo_double_i x);
pseudo_double_i pdi_atan2_rev(pseudo_double_i y, pseudo_double_i x);
pseudo_double_i pdi_sin(pseudo_double_i x);
pseudo_double_i pdi_cos(pseudo_double_i x);
pseudo_double_i pdi_atan2(pseudo_double_i y, pseudo_double_i x);

// x is a 2.62 unsigned fixed in the range (1,4)
// result is 1.63 unsigned fixed in the range (0.5,1)
uint64_t inv_sqrt64_fixed(uint64_t x);

// x is a 0.64 unsigned fixed in the range [0,1)
// result is 2.62 unsigned fixed in the range [1,2)
uint64_t exp2_64_fixed(uint64_t x);

// x is a 1.63 unsigned fixed in the range [0,1)
// calculate ln2(x+1)
// result is 1.63 unsigned fixed in the range [0,1)
uint64_t log2_64_fixed(uint64_t x);

// x is a 2.62 unsigned fixed in the range [0,1]
// calculate sin_rev(x)
// result is 2.62 unsigned fixed in the range [0,1]
uint64_t sin_rev_64_fixed(uint64_t x);

// x is a 2.62 unsigned fixed in the range [0,1]
// calculate atan_rev(x)
// result is 2.62 unsigned fixed in the range [0,1]
uint64_t atan_rev_64_fixed(uint64_t x);

// useful to expose this
inline uint64_t multu64hi(uint64_t x,uint64_t y) {return (((unsigned __int128)x)*y)>>64;}

void debug_pdi_output(pseudo_double_i d);

/* =========================================================================================================================================================================
 * below this line is the C wrapper
 * =========================================================================================================================================================================
 */

struct pseudo_double {
	pseudo_double_i val;
};

inline pseudo_double create_pseudo_double_from_internal(pseudo_double_i x) {
	pseudo_double ret;
	ret.val=x;
	return ret;
}

inline pseudo_double pd_neg(pseudo_double x) {return create_pseudo_double_from_internal(pdi_neg(x.val));}
inline pseudo_double pd_abs(pseudo_double x) {return create_pseudo_double_from_internal(pdi_abs(x.val));}
inline int pd_gt(pseudo_double x, pseudo_double y) {return pdi_gt(x.val,y.val);}
inline int pd_gte(pseudo_double x, pseudo_double y) {return pdi_gte(x.val,y.val);}
inline pseudo_double pd_sub(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_sub(x.val,y.val));}
inline pseudo_double pd_add(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_add(x.val,y.val));}
inline pseudo_double pd_mult(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_mult(x.val,y.val));}
inline pseudo_double pd_div(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_div(x.val,y.val));}
inline pseudo_double pd_ldexp(pseudo_double x, int y) {return create_pseudo_double_from_internal(pdi_ldexp(x.val,y));}
inline pseudo_double double_to_pd(double d) {return create_pseudo_double_from_internal(double_to_pdi(d));}
inline pseudo_double int64_to_pd(int64_t d) {return create_pseudo_double_from_internal(int64_to_pdi(d));}
inline pseudo_double uint64_to_pd(uint64_t d) {return create_pseudo_double_from_internal(uint64_to_pdi(d));}
inline double pd_to_double(pseudo_double d) {return pdi_to_double(d.val);}
inline int64_t pd_to_int64(pseudo_double d) {return pdi_to_int64(d.val);}
inline uint64_t pd_to_uint64(pseudo_double d) {return pdi_to_uint64(d.val);}
inline pseudo_double int64fixed10_to_pd(int64_t d, int32_t e) {return create_pseudo_double_from_internal(int64fixed10_to_pdi(d,e));}
inline pseudo_double int64fixed2_to_pd(int64_t d, int32_t e) {return create_pseudo_double_from_internal(int64fixed2_to_pdi(d,e));}
inline int64_t pd_to_int64fixed2(pseudo_double d, int32_t e) {return pdi_to_int64fixed2(d.val,e);}
inline pseudo_double pd_floor(pseudo_double x) {return create_pseudo_double_from_internal(pdi_floor(x.val));}
inline pseudo_double pd_ceil(pseudo_double x) {return create_pseudo_double_from_internal(pdi_ceil(x.val));}
inline pseudo_double pd_round(pseudo_double x) {return create_pseudo_double_from_internal(pdi_round(x.val));}
inline pseudo_double pd_sqrt(pseudo_double x) {return create_pseudo_double_from_internal(pdi_sqrt(x.val));}
inline pseudo_double pd_inv_sqrt(pseudo_double x) {return create_pseudo_double_from_internal(pdi_inv_sqrt(x.val));}
inline pseudo_double pd_exp2(pseudo_double x) {return create_pseudo_double_from_internal(pdi_exp2(x.val));}
inline pseudo_double pd_exp(pseudo_double x) {return create_pseudo_double_from_internal(pdi_exp(x.val));}
inline pseudo_double pd_log2(pseudo_double x) {return create_pseudo_double_from_internal(pdi_log2(x.val));}
inline pseudo_double pd_log(pseudo_double x) {return create_pseudo_double_from_internal(pdi_log(x.val));}
inline pseudo_double pd_log10(pseudo_double x) {return create_pseudo_double_from_internal(pdi_log10(x.val));}
inline pseudo_double pd_pow(pseudo_double x, pseudo_double y) {return create_pseudo_double_from_internal(pdi_pow(x.val,y.val));}
inline pseudo_double pd_sin_rev(pseudo_double x) {return create_pseudo_double_from_internal(pdi_neg(x.val));}
inline pseudo_double pd_cos_rev(pseudo_double x) {return create_pseudo_double_from_internal(pdi_cos_rev(x.val));}
inline pseudo_double pd_atan2_rev(pseudo_double y, pseudo_double x) {return create_pseudo_double_from_internal(pdi_atan2_rev(x.val,y.val));}
inline pseudo_double pd_sin(pseudo_double x) {return create_pseudo_double_from_internal(pdi_sin(x.val));}
inline pseudo_double pd_cos(pseudo_double x) {return create_pseudo_double_from_internal(pdi_cos(x.val));}
inline pseudo_double pd_atan2(pseudo_double y, pseudo_double x) {return create_pseudo_double_from_internal(pdi_atan2(x.val,y.val));}
inline void debug_pd_output(pseudo_double d) {debug_pdi_output(d.val);}

#endif // PSEUDO_DOUBLE_H
