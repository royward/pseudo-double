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
#include <math.h>
#include <stdio.h>

pseudo_double double_to_pd(double d) {
	union {
		uint64_t i;
		double d;
	} v;
	if(d==0) {
		return 0;
	}
	v.d=d;
	uint64_t i=v.i;
	int negative=(((int64_t)i)<0);
	int32_t raw_exponent=((i>>52)&0x7FF);
	signed_pd_internal exponent=raw_exponent+PSEUDO_DOUBLE_EXP_BIAS-0x3FF+2;
	int64_t old_mantissa=(i&0xFFFFFFFFFFFFFULL);
	signed_pd_internal mantissa=old_mantissa+0x10000000000000ULL; // add in the implied bit
	if(negative) {
		if(old_mantissa==0) {
			if(exponent<1) {
				return 0;
			}
#if PF_ERROR_CHECK
			if(exponent>(signed_pd_internal)(EXP_MASK+1)) {
				PF_DO_ERROR_OVERFLOW;
			}
#endif
			return (1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-1))+exponent-1;
		}
	}
#if PF_ERROR_CHECK
	if(exponent<0) {
		return 0;
	}
	if(exponent>(signed_pd_internal)EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
#endif
	mantissa=shift_left_signed(mantissa,PSEUDO_DOUBLE_TOTAL_BITS-54);
	//mantissa=(mantissa+PSEUDO_DOUBLE_HALF_ULP)&~PSEUDO_DOUBLE_HALF_ULP;
	if(negative) {
		return -(mantissa&EXP_MASK_INV)+exponent;
	} else {
		return (mantissa&EXP_MASK_INV)+exponent;
	}
}

pseudo_double int64fixed2_to_pd(int64_t d, int32_t e) {
	if(d==0) {
		return 0;
	}
	int negative=(d<0);
	int lead_bits=clz(negative?~d:d);
	return ((shift_left_signed(d,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits+e;
}

pseudo_double int64_to_pd(int64_t d) {
	if(d==0) {
		return 0;
	}
	int negative=(d<0);
	int lead_bits=clz(negative?~d:d);
	return ((shift_left_signed(d,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits;
}

pseudo_double uint64_to_pd(uint64_t d) {
	if(d==0) {
		return 0;
	}
	int lead_bits=clz(d);
	return ((shift_left_signed(d,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits;
}

double pd_to_double(pseudo_double x) {
	union {
		uint64_t i;
		double d;
	} v;
	if(x==0) {
		return 0.0;
	}
	signed_pd_internal vx=((signed_pd_internal)x)&EXP_MASK_INV;
	uint64_t sgn=0;
	int32_t exponent=(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS+0x3FF-2;
	if(vx<0) {
		sgn=0x8000000000000000ULL;
		if(((uint64_t)vx)==(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-1))) {
			if(exponent<-1) {
				return 0;
			}
			if(exponent>=0x7FF) {
				return NAN;
			}
			v.i=(((uint64_t)exponent+1)<<52)+sgn;
			return v.d;
		}
		vx=-vx;
	}
	if(exponent<0) {
		return 0;
	}
	if(exponent>0x7FF) {
		return NAN;
	}
	v.i=((((unsigned_pd_internal)vx)<<2)>>12)+(((uint64_t)exponent)<<52)+sgn;
	return v.d;
}

int64_t pd_to_int64(pseudo_double x) {
	if(x==0) {
		return 0;
	}
	signed_pd_internal vx=((signed_pd_internal)x)&EXP_MASK_INV;
	int32_t exponent=(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
#if PF_ERROR_CHECK
	if(exponent>PSEUDO_DOUBLE_TOTAL_BITS) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(PSEUDO_DOUBLE_TOTAL_BITS-exponent>=64) {
		return 0;
	}
#endif
	int64_t ret=vx>>(PSEUDO_DOUBLE_TOTAL_BITS-exponent);
	return ret;
}

int64_t pd_to_int64fixed2(pseudo_double x, int32_t e) {
	if(x==0) {
		return 0;
	}
	signed_pd_internal vx=((signed_pd_internal)x)&EXP_MASK_INV;
	int32_t exponent=(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS-e;
#if PF_ERROR_CHECK
	if(exponent>PSEUDO_DOUBLE_TOTAL_BITS) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(PSEUDO_DOUBLE_TOTAL_BITS-exponent>=64) {
		PF_DO_ERROR_UNDERFLOW;
	}
#endif
	int64_t ret=vx>>(PSEUDO_DOUBLE_TOTAL_BITS-exponent);
	return ret;
}

uint64_t pd_to_uint64(pseudo_double x) {
	if(((signed_pd_internal)x)<0) {
		PF_DO_ERROR_RANGE;
	}
	if(x==0) {
		return 0;
	}
	unsigned_pd_internal vx=((unsigned_pd_internal)x)&EXP_MASK_INV;
	int32_t exponent=(x&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	if(exponent==PSEUDO_DOUBLE_TOTAL_BITS+1) {
		return vx<<1;
	}
#if PF_ERROR_CHECK
	if(exponent>PSEUDO_DOUBLE_TOTAL_BITS) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(PSEUDO_DOUBLE_TOTAL_BITS-exponent>=64) {
		PF_DO_ERROR_UNDERFLOW;
	}
#endif
	uint64_t ret=vx>>(PSEUDO_DOUBLE_TOTAL_BITS-exponent);
	return ret;
}

// x is a 2.62 unsigned fixed in the range (1,4)
// result is 1.63 unsigned fixed in the range (0.5,1)
uint64_t inv_sqrt64_fixed(uint64_t x) {
	// start with a linear interpolation correct at the endpoints
	// 7/6 - 1/6 x, so 1->1, 4->0.5
	uint64_t y=3074457345618258602ULL-multu64hi(x,12297829382473034410ULL);
	// now do some Newton-Raphson
	// y=y*(3/2-1/2*x*y*y)
	// Maximum error for #iterations:
	// 0	~0.2
	// 1	~0.06
	// 2	~0.005
	// 3	~3e-5
	// 4	~1e-9
	// 5	~1e-18 (about 60 bits - limit of the algorithm)
	y=multu64hi(y,0xC000000000000000ULL-((multu64hi(multu64hi(y,y),x))))<<1;
	y=multu64hi(y,0xC000000000000000ULL-((multu64hi(multu64hi(y,y),x))))<<1;
	y=multu64hi(y,0xC000000000000000ULL-((multu64hi(multu64hi(y,y),x))))<<1;
	y=multu64hi(y,0xC000000000000000ULL-((multu64hi(multu64hi(y,y),x))))<<1;
	y=multu64hi(y,0xC000000000000000ULL-((multu64hi(multu64hi(y,y),x)))); // dont shift left on the last one
	return y;
}

/****************************************************************************************************
 * Uncomment this if we need a 32 bit version. Similiar to the 64 bit version, but has not been tested.
uint32_t multu32hi(uint32_t x,uint32_t y) {
	return (((unsigned uint64_t)x)*y)>>32;
}

// x is a 2.30 unsigned fixed in the range (1,4)
// result is 1.31 unsigned fixed in the range (0.5,1)
uint32_t inv_sqrt32_internal(uint32_t x) {
	// start with a linear interpolation correct at the endpoints
	// 7/6 - 1/6 x, so 1->1, 4->0.5
	uint32_t y=715827882U-multu32hi(x,2863311530U);
	// now do some Newton-Raphson
	// y=y*(3/2-1/2*x*y*y)
	// Maximum error for #iterations:
	// 0	~0.2
	// 1	~0.06
	// 2	~0.005
	// 3	~3e-5
	// 4	~1e-9 (about 30 bits - limit of the algorithm)
	y=multu32hi(y,0xC0000000U-((multu32hi(multu32hi(y,y),x))))<<1;
	y=multu32hi(y,0xC0000000U-((multu32hi(multu32hi(y,y),x))))<<1;
	y=multu32hi(y,0xC0000000U-((multu32hi(multu32hi(y,y),x))))<<1;
	y=multu32hi(y,0xC0000000U-((multu32hi(multu32hi(y,y),x)))); // dont shift left on the last one
	return y;
}
****************************************************************************************************/

pseudo_double pd_inv_sqrt(pseudo_double x) {
#if PF_ERROR_CHECK
	if(((signed_pd_internal)x)<=0) {
		PF_DO_ERROR_RANGE;
	}
#endif
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
	return (inv_sqrt64_fixed(mantissa)&EXP_MASK_INV)+3*(PSEUDO_DOUBLE_EXP_BIAS>>1)+2-(exponent>>1);;
}

pseudo_double pd_sqrt(pseudo_double x) {
#if PF_ERROR_CHECK
	if(((signed_pd_internal)x)<0) {
		PF_DO_ERROR_RANGE;
	}
#endif
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
	uint64_t y=multu64hi((inv_sqrt64_fixed(mantissa>>(64-PSEUDO_DOUBLE_TOTAL_BITS))<<(64-PSEUDO_DOUBLE_TOTAL_BITS)),mantissa)<<1;
	//printf("%lx * %lx = %lx\n",inv_sqrt_internal(mantissa),mantissa,y);
	return (y&EXP_MASK_INV)+(PSEUDO_DOUBLE_EXP_BIAS>>1)+1+(exponent>>1);
}

// ./lolremez --long-double -d 10 -r "0:1" "exp2(x)"
// long double f(long double x) {
//     long double u = 1.0006697217452461573e-8l;
//     u = u * x + 9.4339203623592751071e-8l;
//     u = u * x + 1.3318412101751585619e-6l;
//     u = u * x + 1.5243996624121435679e-5l;
//     u = u * x + 1.5404004101871478009e-4l;
//     u = u * x + 1.3333541706066723834e-3l;
//     u = u * x + 9.6181294626261305208e-3l;
//     u = u * x + 5.5504108620159730484e-2l;
//     u = u * x + 2.4022650696198007158e-1l;
//     u = u * x + 6.9314718055987294141e-1l;
//     return u * x + 1.0000000000000003006l;
// }

// x is a 0.64 unsigned fixed in the range [0,1)
// result is 2.62 unsigned fixed in the range [1,2)
uint64_t exp2_64_fixed(uint64_t x) {
	uint64_t u=184590982593ULL;
	u=multu64hi(u,x)+1740251145362ULL;
	u=multu64hi(u,x)+24568133950921ULL;
	u=multu64hi(u,x)+281202104385660ULL;
	u=multu64hi(u,x)+2841537213775953ULL;
	u=multu64hi(u,x)+24596043144794548ULL;
	u=multu64hi(u,x)+177423172664869807ULL;
	u=multu64hi(u,x)+1023870086755462747ULL;
	u=multu64hi(u,x)+4431396893648852228ULL;
	u=multu64hi(u,x)+(12786308645201320706ULL+0x2B5B);
	return (multu64hi(u,x)>>2)+0x4000000000000000ULL;
}

int64_t mults64hir1(int64_t x,int64_t y) {
	return ((((__int128)x)*y)>>64)<<1;
}

int64_t mults64hi(int64_t x,int64_t y) {
	return ((((__int128)x)*y)>>64);
}

// ./lolremez --stats --debug --long-double -d 18 -r "0:1" "log2(x+1)"
// long double f(long double x) {
//     long double u = -9.3911951398832967647e-5l;
//     u = u * x + 9.8625752303230521983e-4l;
//     u = u * x + -4.9037756240706259112e-3l;
//     u = u * x + 1.5466003256993887004e-2l;
//     u = u * x + -3.5114005964259793157e-2l;
//     u = u * x + 6.2097710142540749261e-2l;
//     u = u * x + -9.1018203799855508794e-2l;
//     u = u * x + 1.1695926184602648963e-1l;
//     u = u * x + -1.3875073757264520723e-1l;
//     u = u * x + 1.5861010711744948314e-1l;
//     u = u * x + -1.7993688675087083643e-1l;
//     u = u * x + 2.0602762670500648349e-1l;
//     u = u * x + -2.4043973440620491453e-1l;
//     u = u * x + 2.8853812976087160084e-1l;
//     u = u * x + -3.6067370565679917042e-1l;
//     u = u * x + 4.8089834488826342563e-1l;
//     u = u * x + -7.2134752040270974771e-1l;
//     u = u * x + 1.4426950408886293243l;
//     return u * x;
// }

// x is a 1.63 unsigned fixed in the range [0,1)
// calculate ln2(x+1)
// result is 1.63 unsigned fixed in the range [0,1)
uint64_t log2_64_fixed(uint64_t x) {
	int64_t u=              -866184866458461LL*256;
	u=mults64hi(u,x)       +9096620059073819LL*128;
	u=mults64hi(u,x)      -45229346966063088LL*64;
	u=mults64hi(u,x)     +142648701962462304LL*32;
	u=mults64hi(u,x)     -323869540712705594LL*16;
	u=mults64hi(u,x)     +572750283281423541LL*8;
	u=mults64hi(u,x)     -839494755772336399LL*4;
	u=mults64hi(u,x)    +1078758785161816410LL*2;
	u=mults64hi(u,x)    -1279749673020511097LL;
	u=mults64hi(u<<2,x) +1462920026749624213LL*2;
	u=mults64hi(u,x)    -1659624849656686669LL;
	u=mults64hi(u<<2,x) +1900269450970511052LL*2;
	u=mults64hi(u,x)    -2217665122870979542LL;
	u=mults64hi(u<<1,x) +2661294517602797903LL;
	u=mults64hi(u<<1,x) -3326627771183711640LL;
	u=mults64hi(u<<1,x) +4435504346812152696LL;
	u=mults64hi(u<<1,x) -6653256548536882955LL;
	u=mults64hi(u,x)+   (6653256548920620560LL+0x1005);
//                       9223372036854775808
	return mults64hi(u,x)<<2;
}

// x^y = e^ln(x)^y = e^(y*ln(x))
pseudo_double pd_pow(pseudo_double x, pseudo_double y) {
#if PF_ERROR_CHECK
	if(((signed_pd_internal)x)<=0) {
		PF_DO_ERROR_RANGE;
	}
#endif
	int64_t exponent=(x&EXP_MASK);
	int64_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS-2;
	uint64_t mantissa=((x&EXP_MASK_INV)<<2)>>1;
	uint64_t log_frac=log2_64_fixed(mantissa);
	int64_t vx;
	int64_t expx;
	if(e==0) {
		if(log_frac==0) {
			return uint64_to_pd(1); // 1^z=1
		}
		int lead_bits=clz(log_frac);
		vx=log_frac<<(lead_bits-1);
		expx=2-lead_bits;
	} else if(e==-1) {
		log_frac+=0x8000000000000000ULL;
		int lead_bits=clz(~log_frac);
		vx=log_frac<<(lead_bits-1);
		expx=2-lead_bits;
	} else {
		int negative=(e<0);
		int lead_bits=clz(negative?~e:e);
		vx=((e<<(PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))+(log_frac>>(64-lead_bits)));
		expx=65-lead_bits;
	}
	int32_t expy=(y&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal vy=(signed_pd_internal)(y&EXP_MASK_INV);
	signed_pd_internal vr=(((signed_large_pd_internal)vx)*vy)>>64;
	if(vr==0) {
		// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
		return uint64_to_pd(1); // 2^0=1
	}
	int32_t leading_bits=clz(vr>0?vr:~vr)-1;
	vr<<=leading_bits;
	int32_t er=expx+expy-leading_bits;
	int32_t new_exponent;
	uint64_t fraction;
	if(er<2) {
		if(vr<0) {
			new_exponent=-1;
			if(er==1) {
				fraction=vr<<1;
			} else {
				fraction=vr>>-er;
			}
		} else {
			new_exponent=0;
			fraction=vr;
			if(er==1) {
				fraction<<=1;
			} else {
				fraction>>=-er;
			}
		}
	} else if(er<=PSEUDO_DOUBLE_EXP_BITS) { // max=2^(2^PSEUDO_DOUBLE_EXP_BITS)), log2(max)=2^PSEUDO_DOUBLE_EXP_BITS
		uint64_t m=(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-er))-1;
		new_exponent=((signed_pd_internal)(vr&~m))>>(PSEUDO_DOUBLE_TOTAL_BITS-er);
		fraction=((signed_pd_internal)(vr&m))<<er;
	} else {
#if PF_ERROR_CHECK
		if(vr<0) {
			PF_DO_ERROR_UNDERFLOW;
		} else {
			PF_DO_ERROR_OVERFLOW;
		}
#endif
	}
	int32_t newe=new_exponent+PSEUDO_DOUBLE_EXP_BIAS+2;
	if(newe<0) { // common to have underflow, so leave this in
		PF_DO_ERROR_UNDERFLOW;
#if PF_ERROR_CHECK
	} else if(newe>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
#endif
	}
	//printf("%16lx:%5.10f:%d:%16lx:%f\n",x,pd_to_double(x),new_exponent,fraction,ldexp(fraction,-64));
	return newe+((exp2_64_fixed(fraction<<(64-PSEUDO_DOUBLE_TOTAL_BITS))<<(64-PSEUDO_DOUBLE_TOTAL_BITS))&EXP_MASK_INV);

}

pseudo_double pd_log2(pseudo_double x) {
#if PF_ERROR_CHECK
	if(((signed_pd_internal)x)<=0) {
		PF_DO_ERROR_RANGE;
	}
#endif
	int64_t exponent=(x&EXP_MASK);
	int64_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS-2;
	uint64_t mantissa=((x&EXP_MASK_INV)<<2)>>1;
	uint64_t log_frac=log2_64_fixed(mantissa);
	if(e==0) {
		if(log_frac==0) {
			return 0;
		}
		int lead_bits=clz(log_frac);
		return ((log_frac<<(lead_bits-1))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+2-lead_bits;
	} else if(e==-1) {
		log_frac+=0x8000000000000000ULL;
		int lead_bits=clz(~log_frac);
		return ((log_frac<<(lead_bits-1))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+2-lead_bits;
	}
	int negative=(e<0);
	int lead_bits=clz(negative?~e:e);
	return (((e<<(PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))+(log_frac>>(64-lead_bits)))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits;
// 	printf("%16lx:%5.10f:%5.10f:%d:%16lx:%5.10f:%5.10f\n",x,pd_to_double(x),log2(pd_to_double(x)),e,mantissa,ldexp(1.0+ldexp(mantissa,-63),e),ldexp(log2_64_fixed(mantissa),-63));
// 	return 0;
}

// e^x=(2^log2(e))^x=2^(log2(e)*x)
pseudo_double pd_exp2(pseudo_double x) {
	if(x==0) {
		return uint64_to_pd(1);
	}
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	int32_t new_exponent;
	uint64_t fraction;
	if(e<2) {
		if(((signed_pd_internal)x)<0) {
			new_exponent=-1;
			if(e==1) {
				fraction=(x&EXP_MASK_INV)<<1;
			} else {
				fraction=((signed_pd_internal)(x&EXP_MASK_INV))>>-e;
			}
		} else {
			new_exponent=0;
			fraction=(x&EXP_MASK_INV);
			if(e==1) {
				fraction<<=1;
			} else {
				fraction>>=-e;
			}
		}
	} else if(e<=PSEUDO_DOUBLE_EXP_BITS) { // max=2^(2^PSEUDO_DOUBLE_EXP_BITS)), log2(max)=2^PSEUDO_DOUBLE_EXP_BITS
		uint64_t m=(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-e))-1;
		new_exponent=((signed_pd_internal)(x&~m))>>(PSEUDO_DOUBLE_TOTAL_BITS-e);
		fraction=((signed_pd_internal)(x&EXP_MASK_INV&m))<<e;
	} else {
		// common to have underflow, so leave this in
		if(((signed_pd_internal)x)<0) {
			PF_DO_ERROR_UNDERFLOW;
#if PF_ERROR_CHECK
		} else {
			PF_DO_ERROR_OVERFLOW;
#endif
		}
	}
	int32_t newe=new_exponent+PSEUDO_DOUBLE_EXP_BIAS+2;
#if PF_ERROR_CHECK
	if(newe<0) {
		PF_DO_ERROR_UNDERFLOW;
	} else if(newe>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
#endif
	//printf("%16lx:%5.10f:%d:%16lx:%f\n",x,pd_to_double(x),new_exponent,fraction,ldexp(fraction,-64));
	return newe+((exp2_64_fixed(fraction<<(64-PSEUDO_DOUBLE_TOTAL_BITS))<<(64-PSEUDO_DOUBLE_TOTAL_BITS))&EXP_MASK_INV);
}

pseudo_double pd_exp(pseudo_double x) {
	const pseudo_double log_2_e=int64fixed10_to_pd(1442695040888963407,-18);
	return pd_exp2(pd_mult(x,log_2_e));
}

static const pseudo_double pd_inv_log_2_e =int64fixed10_to_pd(6931471805599453094,-19);
static const pseudo_double pd_inv_log_2_10=int64fixed10_to_pd(3010299956639811952,-19);
static const pseudo_double pd_1_div_tau=pd_div(uint64_to_pd(1),int64fixed10_to_pd(6283185307179586477,-18));
static const pseudo_double pd_tau=int64fixed10_to_pd(6283185307179586477,-18);

pseudo_double pd_log(pseudo_double x) {
	return pd_mult(pd_log2(x),pd_inv_log_2_e);
}

pseudo_double pd_log10(pseudo_double x) {
	return pd_mult(pd_log2(x),pd_inv_log_2_10);
}

// ./lolremez --stats --debug --long-double -d 15 -r "-1:1" "sin(x*pi/2)"
// long double f(long double x) {
//     x2 = x * x;
//     long double u = -6.4348418488777433292e-10;
//     u = u * x2 + 5.6877693955975479363e-8;
//     u = u * x2 + -3.5988024751262635033e-6;
//     u = u * x2 + 1.6044116327434520208e-4;
//     u = u * x2 + -4.6817541288738726714e-3;
//     u = u * x2 + 7.9692626245142861354e-2;
//     u = u * x2 + -6.4596409750617316605e-1;
//     u = u * x2 + 1.5707963267948950979l;
//     return u * x;
// }

// x is a 2.62 unsigned fixed in the range [0,1]
// calculate sin_rev_64_fixed(x)
// result is 2.62 unsigned fixed in the range [0,1]
uint64_t sin_rev_64_fixed(uint64_t x) {
    int64_t x2=mults64hi(x,x)<<2;
    int64_t u=                     -2967547018LL;
    u=mults64hi(u<<2,x2)         +262302065977LL;
    u=mults64hi(u<<2,x2)       -16596547057622LL;
    u=mults64hi(u<<2,x2)      +739904269452523LL;
    u=mults64hi(u<<2,x2)    -21590780057842334LL;
    u=mults64hi(u<<2,x2)   +367517370226484839LL;
    u=mults64hi(u<<2,x2)  -2978983596875284700LL;
    u=(mults64hi(u,x2)<<2)+(7244019458077115826LL+415);
    return mults64hi(u,x)<<2;
}

pseudo_double pd_sin_rev(pseudo_double x) {
	if(x==0) {
		return uint64_to_pd(0);
	}
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	uint64_t fraction;
	if(e<2) {
		if(((signed_pd_internal)x)<0) {
			if(e==1) {
				fraction=(x&EXP_MASK_INV)<<1;
			} else {
				fraction=((signed_pd_internal)(x&EXP_MASK_INV))>>-e;
			}
		} else {
			fraction=(x&EXP_MASK_INV);
			if(e==1) {
				fraction<<=1;
			} else {
				fraction>>=-e;
			}
		}
	} else {
		uint64_t m=(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-e))-1;
		fraction=((signed_pd_internal)(x&EXP_MASK_INV&m))<<e;
	}
	int negative=((signed_pd_internal)fraction<0);
	if(negative) {
		fraction=-(signed_pd_internal)fraction;
	}
	if(fraction>>62) {
		fraction=0x8000000000000000ULL-fraction;
	}
	int64_t d=sin_rev_64_fixed(fraction);
	if(d==0) {
		return 0;
	}
	if(negative) {
		d=-d;
	}
	int lead_bits=clz(negative?~d:d);
	return (shift_left_signed(d,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65)&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+3-lead_bits;
}

pseudo_double pd_cos_rev(pseudo_double x) {
	if(x==0) {
		return uint64_to_pd(1);
	}
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	uint64_t fraction;
	if(e<2) {
		if(((signed_pd_internal)x)<0) {
			if(e==1) {
				fraction=(x&EXP_MASK_INV)<<1;
			} else {
				fraction=((signed_pd_internal)(x&EXP_MASK_INV))>>-e;
			}
		} else {
			fraction=(x&EXP_MASK_INV);
			if(e==1) {
				fraction<<=1;
			} else {
				fraction>>=-e;
			}
		}
	} else {
		uint64_t m=(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-e))-1;
		fraction=((signed_pd_internal)(x&EXP_MASK_INV&m))<<e;
	}
	fraction+=0x4000000000000000ULL; // _only line that is different between sin and cos
	int negative=((signed_pd_internal)fraction<0);
	if(negative) {
		fraction=-(signed_pd_internal)fraction;
	}
	if(fraction>>62) {
		fraction=0x8000000000000000ULL-fraction;
	}
	int64_t d=sin_rev_64_fixed(fraction);
	if(d==0) {
		return 0;
	}
	if(negative) {
		d=-d;
	}
	int lead_bits=clz(negative?~d:d);
	return (shift_left_signed(d,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65)&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+3-lead_bits;
}

pseudo_double pd_sin(pseudo_double x) {
	return pd_sin_rev(pd_mult(x,pd_1_div_tau));
}

pseudo_double pd_cos(pseudo_double x) {
	return pd_cos_rev(pd_mult(x,pd_1_div_tau));
}

// ./lolremez --stats --debug --long-double -d 35 -r "-1:1" "atan(x)*4/pi"
// long double f(long double x) {
//     x2 = x * x;
//     u = -5.1448538374132849052e-5;
//     u = u * x2 + 5.275833184610939186e-4;
//     u = u * x2 + -2.5593891803515647842e-3;
//     u = u * x2 + 7.8734531693353488079e-3;
//     u = u * x2 + -1.744662472591652788e-2;
//     u = u * x2 + 3.0173042794079607078e-2;
//     u = u * x2 + -4.3417549325196383149e-2;
//     u = u * x2 + 5.5090280539729516217e-2;
//     u = u * x2 + -6.4974658268428946492e-2;
//     u = u * x2 + 7.4312120990007820621e-2;
//     u = u * x2 + -8.475459097263842221e-2;
//     u = u * x2 + 9.7920586959909512452e-2;
//     u = u * x2 + -1.1574658642078503641e-1;
//     u = u * x2 + 1.4147086102567709565e-1;
//     u = u * x2 + -1.8189135316383373911e-1;
//     u = u * x2 + 2.5464790863720206134e-1;
//     u = u * x2 + -4.2441318157402211219e-1;
//     u = u * x2 + 1.2732395447351443229;
//     return u * x;
// }

// x is a 2.62 unsigned fixed in the range [0,1]
// calculate atan_rev_64_fixed(x)
// result is 2.62 unsigned fixed in the range [0,1]
uint64_t atan_rev_64_fixed(uint64_t x) {
	int64_t x2=mults64hi(x,x)<<2;
	int64_t u=                 -237264505088513LL;
	u=mults64hi(u<<2,x2)      +2433048613302551LL;
	u=mults64hi(u<<2,x2)     -11803099298741644LL;
	u=mults64hi(u<<2,x2)     +36309893897766633LL;
	u=mults64hi(u<<2,x2)     -80458355317258810LL;
	u=mults64hi(u<<2,x2)    +139148599586868171LL;
	u=mults64hi(u<<2,x2)    -200228105177389631LL;
	u=mults64hi(u<<2,x2)    +254059076516313023LL;
	u=mults64hi(u<<2,x2)    -299642723088611246LL;
	u=mults64hi(u<<2,x2)    +342704169369303486LL;
	u=mults64hi(u<<2,x2)    -390861562186048719LL;
	u=mults64hi(u<<2,x2)    +451579001799217900LL;
	u=mults64hi(u<<2,x2)    -533786914277431708LL;
	u=mults64hi(u<<2,x2)    +652419191806999136LL;
	u=mults64hi(u<<2,x2)    -838825810258490282LL;
	u=mults64hi(u<<2,x2)   +1174356199883959617LL;
	u=mults64hi(u<<2,x2)   -1957260335501202067LL;
	u=(mults64hi(u,x2)<<2)+(5871781006563917768LL+2243);
    return mults64hi(u,x)<<2;
}

pseudo_double pd_atan2_rev(pseudo_double y, pseudo_double x) {
	int negative=0; // boolean
	uint64_t add_const;
	if(y==0) {
		if(((signed_pd_internal)x)>=0) {
			return uint64_to_pd(0);
		} else {
			return int64fixed2_to_pd(1,-1); // 1/2
		}
	} else if(((signed_pd_internal)y)>0) {
		if(x==0) {
			return int64fixed2_to_pd(1,-2); // 1/4
		} else if (((signed_pd_internal)x)>0) {
			if(pd_gte(x,y)) {
				// q1
				add_const=0;
			} else {
				// q2
				{ pseudo_double t=y; y=x; x=t;}
				add_const=0x4000000000000000ULL;
				negative=1; // boolean
			}
		} else { // x<0
			x=pd_neg(x);
			if(pd_gte(x,y)) {
				// q4
				add_const=0x8000000000000000ULL;
				negative=1; // boolean
			} else {
				// q3
				{ pseudo_double t=y; y=x; x=t;}
				add_const=0x4000000000000000ULL;
			}
		}
	} else { // y<0
		y=pd_neg(y);
		if(x==0) {
			return int64fixed2_to_pd(3,-2); // 3/2
		} else if (((signed_pd_internal)x)>0) {
			if(pd_gte(x,y)) {
				// q8
				add_const=0; // boolean
				negative=1;
			} else {
				// q7
				{ pseudo_double t=y; y=x; x=t;}
				add_const=0xC000000000000000ULL;
			}
		} else { // x<0
			x=pd_neg(x);
			if(pd_gte(x,y)) {
				// q5
				add_const=0x8000000000000000ULL;
			} else {
				// q6
				{ pseudo_double t=y; y=x; x=t;}
				add_const=0xC000000000000000ULL;
				negative=1; // boolean
			}
		}
	}
	int32_t expx=x&EXP_MASK;
	int32_t expy=y&EXP_MASK;
	signed_pd_internal vx=(signed_pd_internal)(x&EXP_MASK_INV);
	signed_pd_internal vy=(signed_pd_internal)(y&EXP_MASK_INV);
	uint64_t ratio;
	if(x==y) {
		ratio=uint64_to_pd(1);
	} else if(x==pd_neg(y)) {
		ratio=int64_to_pd(-1);
	} else {
		signed_pd_internal vr=(((((signed_large_pd_internal)vy)>>2)<<64)/vx);
		if(vr==0) {
			ratio=0;
		} else {
			int32_t leading_bits=clz(vr)-1;
			vr<<=leading_bits;
			int32_t new_exponent=expy-expx-leading_bits;
			if(new_exponent<-63) {
				ratio=0;
			} else {
				ratio=vr>>(-new_exponent);
			}
		}
	}
	int64_t d=atan_rev_64_fixed(ratio)>>1;
	d=add_const+(negative?-d:d);
	if(d==0) {
		return 0;
	}
	int negatived=(d<0);
	int lead_bits=clz(negatived?~d:d);
	return ((shift_left_signed(d,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+1-lead_bits;
}

pseudo_double pd_atan2(pseudo_double y, pseudo_double x) {
	const pseudo_double pd_tau=int64fixed10_to_pd(6283185307179586477,-18);
	return pd_mult(pd_atan2_rev(y,x),pd_tau);
}

pseudo_double pd_floor(pseudo_double x) {
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	if(e<2) {
		return (((signed_pd_internal)x)<0)?int64_to_pd(-1):0;
	}
	if(e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS) {
		return x;
	}
	uint64_t m=(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-e))-1;
	return (x&~m)+exponent;
}

pseudo_double pd_ceil(pseudo_double x) {
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal mantissa=x&EXP_MASK_INV;
	if(e<2) {
		if(e==1 && (mantissa<<1)==0) { // special test for ceil(-1)=-1
			return x;
		}
		return (((signed_pd_internal)x)>0)?int64_to_pd(1):0;
	}
	if(e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS) {
		return x;
	}
	uint64_t m=(1ULL<<(PSEUDO_DOUBLE_TOTAL_BITS-e-1))-1;
	signed_pd_internal vr=((mantissa>>1)+m)&~m;
	int32_t leading_bits=clz(vr>0?vr:~vr)-1;
	vr<<=leading_bits;
	int32_t new_exponent=exponent+1-leading_bits;
#if PF_ERROR_CHECK
	if(new_exponent>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(new_exponent<0) {
		PF_DO_ERROR_UNDERFLOW;
	}
#endif
	return (pseudo_double)(vr+new_exponent);
}

pseudo_double pd_round(pseudo_double x) {
	int32_t exponent=(x&EXP_MASK);
	int32_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS;
	signed_pd_internal mantissa=x&EXP_MASK_INV;
	if(e<1) {
		if(e==0 && (mantissa<<1)==0) { // special test for round(-0.5)=-1
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
	int32_t leading_bits=clz(vr>0?vr:~vr)-1;
	vr<<=leading_bits;
	int32_t new_exponent=exponent+1-leading_bits;
#if PF_ERROR_CHECK
	if(new_exponent>EXP_MASK) {
		PF_DO_ERROR_OVERFLOW;
	}
	if(new_exponent<0) {
		PF_DO_ERROR_UNDERFLOW;
	}
#endif
	return (pseudo_double)(vr+new_exponent);
}

