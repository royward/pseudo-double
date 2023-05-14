#ifndef PSEUDO_FLOAT_H
#define PSEUDO_FLOAT_H
#include <stdint.h>

// Can set any of these or just use the defaults
//#define PSEUDO_FLOAT_TOTAL_BITS 64
//#define PSEUDO_FLOAT_EXP_BITS 16
// This will default to 2^(PSEUDO_FLOAT_EXP_BITS-1)
//#define PSEUDO_FLOAT_EXP_BIAS 32768
// This will default to return PF_NAN
//#define PF_DO_ERROR_OVERFLOW return PF_NAN
//#define PF_DO_ERROR_UNDERFLOW return 0
//#define PF_DO_ERROR_RANGE return PF_NAN

// Mantissa
// [.1000 .. .1011] = [-0.5 .. -0.25)
// [.0100 .. .0111] = [0.25 .. 0.5)

#ifndef PSEUDO_FLOAT_TOTAL_BITS
#define PSEUDO_FLOAT_TOTAL_BITS 64
#endif

#ifndef PSEUDO_FLOAT_EXP_BITS
// Recommend 8,16 or 32 for this. Other values will work but be less efficient
#define PSEUDO_FLOAT_EXP_BITS 16
#endif

#define EXP_MASK ((1LL<<PSEUDO_FLOAT_EXP_BITS)-1)
#define EXP_MASK_INV (~((1ULL<<PSEUDO_FLOAT_EXP_BITS)-1))
#if PSEUDO_FLOAT_TOTAL_BITS==64
typedef uint64_t pseudo_float;
typedef uint64_t unsigned_pf_internal;
typedef int64_t signed_pf_internal;
typedef __int128 signed_large_pf_internal;
#define clz __builtin_clzll
#elif PSEUDO_FLOAT_TOTAL_BITS==32
typedef uint32_t pseudo_float;
typedef uint32_t unsigned_pf_internal;
typedef int32_t signed_pf_internal;
typedef int64_t signed_large_pf_internal;
#define clz __builtin_clz
#else
#error PSEUDO_FLOAT_TOTAL_BITS must be 32 or 64 bits
#endif
#ifndef PSEUDO_FLOAT_EXP_BIAS
#define PSEUDO_FLOAT_EXP_BIAS (1U<<(PSEUDO_FLOAT_EXP_BITS-1))
#endif
#define PSEUDO_FLOAT_HALF_ULP ((1ULL<<(PSEUDO_FLOAT_EXP_BITS-1))-1)
#define PF_NAN ((pseudo_float)(-1))

#ifndef PF_DO_ERROR_OVERFLOW
#define PF_DO_ERROR_OVERFLOW return PF_NAN
#endif
#ifndef PF_DO_ERROR_UNDERFLOW
#define PF_DO_ERROR_UNDERFLOW return 0
#endif
#ifndef PF_DO_ERROR_RANGE
#define PF_DO_ERROR_RANGE return PF_NAN
#endif

inline signed_pf_internal shift_left_signed(signed_pf_internal x, int shift) {
	if(shift>=0) {
		return x<<shift;
	}
	return x>>-shift;
}

inline unsigned_pf_internal shift_left_unsigned(unsigned_pf_internal x, int shift) {
	if(shift>=0) {
		return x<<shift;
	}
	return x>>-shift;
}

inline pseudo_float pf_neg(pseudo_float x) {
	// check for special cases due to the representation range of two's complement being asymmetric:
	// mantissa of 1000000000... can't be negated directly and will need an increase in the exponent
	// mantissa of 0100000000... has exponent decreased after negation
	// everything else negates as expected
	uint32_t exponent=x&EXP_MASK;
	unsigned_pf_internal mantissa=x&EXP_MASK_INV;
	if((mantissa<<2)==0) {
		// get the high order byte
		uint32_t hi_byte=x>>(PSEUDO_FLOAT_TOTAL_BITS-8);
		if(hi_byte==0x80) {
			if(exponent==EXP_MASK) {
				PF_DO_ERROR_OVERFLOW;
			}
			return (mantissa>>1)+exponent+1;
		}
		if(hi_byte==0x40) {
			if(exponent==0) {
				PF_DO_ERROR_UNDERFLOW;
			}
			return (mantissa<<1)+exponent-1;
		}
	}
	return (-(x&EXP_MASK_INV))+exponent;
}

inline pseudo_float pf_abs(pseudo_float x) {
	if(((signed_pf_internal)x)>=0) {
		return x;
	}
	uint32_t exponent=x&EXP_MASK;
	unsigned_pf_internal mantissa=x&EXP_MASK_INV;
	if((mantissa<<2)==0) {
		// get the high order byte
		uint32_t hi_byte=x>>(PSEUDO_FLOAT_TOTAL_BITS-8);
		if(hi_byte==0x80) {
			if(exponent==EXP_MASK) {
				PF_DO_ERROR_OVERFLOW;
			}
			return (mantissa>>1)+exponent+1;
		}
	}
	return (-(x&EXP_MASK_INV))+exponent;
}

inline int pf_gt(pseudo_float x, pseudo_float y) {
	int neg=((unsigned_pf_internal)y)>>(PSEUDO_FLOAT_TOTAL_BITS-1);
	if((x^y)>>(PSEUDO_FLOAT_TOTAL_BITS-1)) {
		return neg;
	}
	// signs are the same, check exponent
	int expdiff=(signed_pf_internal)((x&EXP_MASK)-(y&EXP_MASK));
	if(expdiff!=0) {
		return (expdiff>0)^neg;
	}
	// exponents are the same, check mantissa
	return (x>y);
}

inline int pf_gte(pseudo_float x, pseudo_float y) {
	int neg=((unsigned_pf_internal)y)>>(PSEUDO_FLOAT_TOTAL_BITS-1);
	if((x^y)>>(PSEUDO_FLOAT_TOTAL_BITS-1)) {
		return neg;
	}
	// signs are the same, check exponent
	int expdiff=(signed_pf_internal)((x&EXP_MASK)-(y&EXP_MASK));
	if(expdiff!=0) {
		return (expdiff>0)^neg;
	}
	// exponents are the same, check mantissa
	return (x>=y);
}

inline pseudo_float pf_sub(pseudo_float x, pseudo_float y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	if(expy-expx>=(PSEUDO_FLOAT_TOTAL_BITS-1)) {
		return pf_neg(y);
	}
	if(expy-expx<=-(PSEUDO_FLOAT_TOTAL_BITS-1)) {
		return x;
	}
	int32_t exp_max=((expy>expx)?expy:expx)+1;
	signed_pf_internal vx=((signed_pf_internal)(x&EXP_MASK_INV))>>(exp_max-expx);
	signed_pf_internal vy=((signed_pf_internal)(y&EXP_MASK_INV))>>(exp_max-expy);
	signed_pf_internal vr=(vx-vy+PSEUDO_FLOAT_HALF_ULP)&~PSEUDO_FLOAT_HALF_ULP;
	if(vr==0) {
		// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
		return (pseudo_float)0;
	}
	int32_t leading_bits=clz(vr>0?vr:~vr)-1;
	if(leading_bits>exp_max) {
		leading_bits=exp_max;
	}
	vr<<=leading_bits;
	int32_t new_exponent=exp_max-leading_bits;
	if(new_exponent>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(new_exponent<0) {
		return 0;
	}
	return (pseudo_float)((vr&EXP_MASK_INV)+new_exponent);
}

inline pseudo_float pf_add(pseudo_float x, pseudo_float y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	if(expy-expx>=(PSEUDO_FLOAT_TOTAL_BITS-1)) {
		return y;
	}
	if(expy-expx<=-(PSEUDO_FLOAT_TOTAL_BITS-1)) {
		return x;
	}
	int32_t exp_max=((expy>expx)?expy:expx)+1;
	signed_pf_internal vx=((signed_pf_internal)(x&EXP_MASK_INV))>>(exp_max-expx);
	signed_pf_internal vy=((signed_pf_internal)(y&EXP_MASK_INV))>>(exp_max-expy);
	signed_pf_internal vr=(vx+vy+PSEUDO_FLOAT_HALF_ULP)&~PSEUDO_FLOAT_HALF_ULP;
	if(vr==0) {
		// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
		return (pseudo_float)0;
	}
	int32_t leading_bits=clz(vr>0?vr:~vr)-1;
	if(leading_bits>exp_max) {
		leading_bits=exp_max;
	}
	vr<<=leading_bits;
	int32_t new_exponent=exp_max-leading_bits;
	if(new_exponent>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(new_exponent<0) {
		return 0;
	}
	return (pseudo_float)((vr&EXP_MASK_INV)+new_exponent);
}

inline pseudo_float pf_mult(pseudo_float x, pseudo_float y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	signed_pf_internal vx=(signed_pf_internal)(x&EXP_MASK_INV);
	signed_pf_internal vy=(signed_pf_internal)(y&EXP_MASK_INV);
	signed_pf_internal vr=(((signed_large_pf_internal)vx)*vy)>>64;
	if(vr==0) {
		// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
		return (pseudo_float)0;
	}
	int32_t leading_bits=clz(vr>0?vr:~vr)-1;
	vr<<=leading_bits;
	int32_t new_exponent=expx+expy-PSEUDO_FLOAT_EXP_BIAS-leading_bits;
	if(new_exponent>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(new_exponent<0) {
		PF_DO_ERROR_UNDERFLOW;
	}
	return (pseudo_float)((vr&EXP_MASK_INV)+new_exponent);
}

inline pseudo_float pf_div(pseudo_float x, pseudo_float y) {
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	signed_pf_internal vx=(signed_pf_internal)(x&EXP_MASK_INV);
	signed_pf_internal vy=(signed_pf_internal)(y&EXP_MASK_INV);
	if(vy==0) {
		PF_DO_ERROR_RANGE;
	}
	signed_pf_internal vr=(((((signed_large_pf_internal)vx)>>2)<<64)/vy);
	if(vr==0) {
		// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
		return (pseudo_float)0;
	}
	int32_t leading_bits=clz(vr>0?vr:~vr)-1;
	vr<<=leading_bits;
	int32_t new_exponent=2+expx-expy+PSEUDO_FLOAT_EXP_BIAS-leading_bits;
	if(new_exponent>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(new_exponent<0) {
		PF_DO_ERROR_UNDERFLOW;
	}
	return (pseudo_float)((vr&EXP_MASK_INV)+new_exponent);
}

inline pseudo_float pf_ldexp(pseudo_float x, int y) {
	int32_t expx=x&EXP_MASK;
	if(expx+y>(int32_t)EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(expx+y<0) {
		PF_DO_ERROR_UNDERFLOW;
	}
	return x+y;
}

pseudo_float double_to_pf(double d);
pseudo_float int64_to_pf(int64_t d);
pseudo_float uint64_to_pf(uint64_t d);
double pf_to_double(pseudo_float d);
int64_t pf_to_int64(pseudo_float d);
uint64_t pf_to_uint64(pseudo_float d);

inline pseudo_float int64fixed10_to_pf(int64_t d, int32_t e) {
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
	int exp=nexp+PSEUDO_FLOAT_EXP_BIAS+65-lead_bits;
	return ((shift_left_signed(d,PSEUDO_FLOAT_TOTAL_BITS+lead_bits-65))&EXP_MASK_INV)+exp;
}

pseudo_float int64fixed2_to_pf(int64_t d, int32_t e);
int64_t pf_to_int64fixed2(pseudo_float d, int32_t e);

pseudo_float pf_floor(pseudo_float x);
pseudo_float pf_ceil(pseudo_float x);
pseudo_float pf_round(pseudo_float x);
pseudo_float pf_sqrt(pseudo_float x);
pseudo_float pf_inv_sqrt(pseudo_float x);
pseudo_float pf_exp2(pseudo_float x);
pseudo_float pf_exp(pseudo_float x);
pseudo_float pf_log2(pseudo_float x);
pseudo_float pf_log(pseudo_float x);
pseudo_float pf_log10(pseudo_float x);
pseudo_float pf_pow(pseudo_float x, pseudo_float y);
pseudo_float pf_sin_rev(pseudo_float x);
pseudo_float pf_cos_rev(pseudo_float x);
pseudo_float pf_atan2_rev(pseudo_float y, pseudo_float x);
pseudo_float pf_sin(pseudo_float x);
pseudo_float pf_cos(pseudo_float x);
pseudo_float pf_atan2(pseudo_float y, pseudo_float x);

// x is a 2.62 unsigned fixed in the range (1,4)
// result is 1.63 unsigned fixed in the range (0.5,1)
uint64_t inv_sqrt64_internal(uint64_t x);

// x is a 0.64 unsigned fixed in the range [0,1)
// result is 2.62 unsigned fixed in the range [1,2)
uint64_t exp2_64_internal(uint64_t x);

// x is a 1.63 unsigned fixed in the range [0,1)
// calculate ln2(x+1)
// result is 1.63 unsigned fixed in the range [0,1)
uint64_t log2_64_internal(uint64_t x);

// x is a 2.62 unsigned fixed in the range [0,1]
// calculate sin_rev_64_internal(x)
// result is 2.62 unsigned fixed in the range [0,1]
uint64_t sin_rev_64_internal(uint64_t x);

// x is a 2.62 unsigned fixed in the range [0,1]
// calculate atan_rev_64_internal(x)
// result is 2.62 unsigned fixed in the range [0,1]
uint64_t atan_rev_64_internal(uint64_t x);

// useful to expose this
inline uint64_t multu64hi(uint64_t x,uint64_t y) {return (((unsigned __int128)x)*y)>>64;}

void debug_pf_output(pseudo_float d);

#endif // PSEUDO_FLOAT_H
