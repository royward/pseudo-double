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

use std::ops::{Add, AddAssign, Sub, SubAssign, Neg, Mul, MulAssign, Div, DivAssign};
use std::cmp::{Eq, Ordering};
use std::convert::From;

#[derive(Debug, Copy, Clone, Eq, PartialEq)]
pub struct PseudoDouble(i64);

const PSEUDO_DOUBLE_TOTAL_BITS: i32 = 64;
const PSEUDO_DOUBLE_EXP_BITS: i32 = 16;
const EXP_MASK: i64 = (1<<PSEUDO_DOUBLE_EXP_BITS)-1;
const EXP_MASK_INV: i64 = !EXP_MASK as i64;
const PSEUDO_DOUBLE_HALF_ULP: i64 = (1<<(PSEUDO_DOUBLE_EXP_BITS-1))-1;
const PSEUDO_DOUBLE_EXP_BIAS: i64 = 1<<(PSEUDO_DOUBLE_EXP_BITS-1);

pub const PD_ZERO: PseudoDouble = PseudoDouble(0);
pub const PD_ONE: PseudoDouble = PseudoDouble::pd_create_mantissa_exp10(1,0);
pub const PD_NEG_ONE: PseudoDouble = PseudoDouble::pd_create_mantissa_exp10(-1,0);
pub const PD_LOG_2_E:      PseudoDouble = PseudoDouble::pd_create_mantissa_exp10(1442695040888963407,-18);
pub const PD_LOG_2_10:     PseudoDouble = PseudoDouble::pd_create_mantissa_exp10(3321928094887362347,-18);
pub const PD_INV_LOG_2_E:  PseudoDouble = PseudoDouble::pd_create_mantissa_exp10(6931471805599453094,-19);
pub const PD_INV_LOG_2_10: PseudoDouble = PseudoDouble::pd_create_mantissa_exp10(3010299956639811952,-19);
pub const PD_TAU:          PseudoDouble = PseudoDouble::pd_create_mantissa_exp10(6283185307179586477,-18);
pub const PD_INV_TAU:      PseudoDouble = PseudoDouble::pd_create_mantissa_exp10(1591549430918953358,-19);

#[inline]
const fn shift_left_signed(x:i64, shift:i32) -> i64 {
	if shift>0 {x<<shift} else {x>>-shift}
}

#[inline]
const fn shift_left_unsigned(x:u64, shift:i32) -> u64 {
	if shift>0 {x<<shift} else {x>>-shift}
}

#[inline]
const fn multu64hi(x:u64, y:u64) -> u64 {
	((x as u128*y as u128)>>64) as u64
}

#[inline]
const fn mults64hi(x:i64, y:i64) -> i64 {
	((x as i128*y as i128)>>64) as i64
}

#[inline]
const fn divs64hi(x:i64, y:i64) -> i64 {
	(((x as i128)<<64)/(y as i128)) as i64
}

// x is a 2.62 unsigned fixed in the range (1,4)
// result is 1.63 unsigned fixed in the range (0.5,1)
const fn inv_sqrt64_fixed(x:u64) -> u64 {
	// start with a linear interpolation correct at the endpoints
	// 7/6 - 1/6 x, so 1->1, 4->0.5
	//let mut y=3074457345618258602u64-multu64hi(x,12297829382473034410u64);
	let mut y=3074457345618258602u64+(-(multu64hi(x,12297829382473034410u64) as i64) as u64);
	// now do some Newton-Raphson
	// y=y*(3/2-1/2*x*y*y)
	// Maximum error for #iterations:
	// 0	~0.2
	// 1	~0.06
	// 2	~0.005
	// 3	~3e-5
	// 4	~1e-9
	// 5	~1e-18 (about 60 bits - limit of the algorithm)
	y=multu64hi(y,0xC000000000000000u64-((multu64hi(multu64hi(y,y),x))))<<1;
	y=multu64hi(y,0xC000000000000000u64-((multu64hi(multu64hi(y,y),x))))<<1;
	y=multu64hi(y,0xC000000000000000u64-((multu64hi(multu64hi(y,y),x))))<<1;
	y=multu64hi(y,0xC000000000000000u64-((multu64hi(multu64hi(y,y),x))))<<1;
	y=multu64hi(y,0xC000000000000000u64-((multu64hi(multu64hi(y,y),x)))); // dont shift left on the last one
	return y;
}

// x is a 0.64 unsigned fixed in the range [0,1)
// result is 2.62 unsigned fixed in the range [1,2)
const fn exp2_64_fixed(x:u64) -> u64 {
	let mut u=184590982593u64;
	u=multu64hi(u,x)+1740251145362u64;
	u=multu64hi(u,x)+24568133950921u64;
	u=multu64hi(u,x)+281202104385660u64;
	u=multu64hi(u,x)+2841537213775953u64;
	u=multu64hi(u,x)+24596043144794548u64;
	u=multu64hi(u,x)+177423172664869807u64;
	u=multu64hi(u,x)+1023870086755462747u64;
	u=multu64hi(u,x)+4431396893648852228u64;
	u=multu64hi(u,x)+(12786308645201320706u64+0x2B5B);
	return (multu64hi(u,x)>>2)+0x4000000000000000u64;
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
const fn log2_64_fixed(xu:u64) -> u64 {
	let x=xu as i64;
	let mut u=              -866184866458461i64*256;
	u=mults64hi(u,x)       +9096620059073819i64*128;
	u=mults64hi(u,x)      -45229346966063088i64*64;
	u=mults64hi(u,x)     +142648701962462304i64*32;
	u=mults64hi(u,x)     -323869540712705594i64*16;
	u=mults64hi(u,x)     +572750283281423541i64*8;
	u=mults64hi(u,x)     -839494755772336399i64*4;
	u=mults64hi(u,x)    +1078758785161816410i64*2;
	u=mults64hi(u,x)    -1279749673020511097i64;
	u=mults64hi(u<<2,x) +1462920026749624213i64*2;
	u=mults64hi(u,x)    -1659624849656686669i64;
	u=mults64hi(u<<2,x) +1900269450970511052i64*2;
	u=mults64hi(u,x)    -2217665122870979542i64;
	u=mults64hi(u<<1,x) +2661294517602797903i64;
	u=mults64hi(u<<1,x) -3326627771183711640i64;
	u=mults64hi(u<<1,x) +4435504346812152696i64;
	u=mults64hi(u<<1,x) -6653256548536882955i64;
	u=mults64hi(u,x)+   (6653256548920620560i64+0x1005);
	return (mults64hi(u,x)<<2) as u64;
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
const fn sin_rev_64_fixed(xu:u64) -> u64 {
	let x=xu as i64;
	let x2=mults64hi(x,x)<<2;
    let mut u=                      -2967547018i64;
    u=mults64hi(u<<2,x2)          +262302065977i64;
    u=mults64hi(u<<2,x2)        -16596547057622i64;
    u=mults64hi(u<<2,x2)       +739904269452523i64;
    u=mults64hi(u<<2,x2)     -21590780057842334i64;
    u=mults64hi(u<<2,x2)    +367517370226484839i64;
    u=mults64hi(u<<2,x2)   -2978983596875284700i64;
    u=(mults64hi(u,x2)<<2)+(7244019458077115826i64+415);
    return (mults64hi(u,x)<<2) as u64;
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
const fn atan_rev_64_fixed(xu:u64) -> u64 {
	let x=xu as i64;
	let x2=mults64hi(x,x)<<2;
	let mut u=                 -237264505088513i64;
	u=mults64hi(u<<2,x2)      +2433048613302551i64;
	u=mults64hi(u<<2,x2)     -11803099298741644i64;
	u=mults64hi(u<<2,x2)     +36309893897766633i64;
	u=mults64hi(u<<2,x2)     -80458355317258810i64;
	u=mults64hi(u<<2,x2)    +139148599586868171i64;
	u=mults64hi(u<<2,x2)    -200228105177389631i64;
	u=mults64hi(u<<2,x2)    +254059076516313023i64;
	u=mults64hi(u<<2,x2)    -299642723088611246i64;
	u=mults64hi(u<<2,x2)    +342704169369303486i64;
	u=mults64hi(u<<2,x2)    -390861562186048719i64;
	u=mults64hi(u<<2,x2)    +451579001799217900i64;
	u=mults64hi(u<<2,x2)    -533786914277431708i64;
	u=mults64hi(u<<2,x2)    +652419191806999136i64;
	u=mults64hi(u<<2,x2)    -838825810258490282i64;
	u=mults64hi(u<<2,x2)   +1174356199883959617i64;
	u=mults64hi(u<<2,x2)   -1957260335501202067i64;
	u=(mults64hi(u,x2)<<2)+(5871781006563917768i64+2243);
    return (mults64hi(u,x)<<2) as u64;
}

impl From<i64> for PseudoDouble {
	fn from(x : i64) -> Self {
		if x==0 {
			return PD_ZERO;
		} else {
			let lead_bits=(if x<0 {!x} else {x}).leading_zeros() as i32;
			return PseudoDouble(((shift_left_signed(x,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits as i64);
		}
    }
}

impl From<u64> for PseudoDouble {
    fn from(x : u64) -> Self {
		if x==0 {
			return PD_ZERO;
		} else {
			let lead_bits=x.leading_zeros() as i32;
			return PseudoDouble(((shift_left_unsigned(x,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65)) as i64&EXP_MASK_INV) +PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits as i64);
		}
    }
}

impl From<PseudoDouble> for i64 {
    fn from(x: PseudoDouble) -> Self {
		if x.0==0 {
			return 0;
		}
		let exponent=((x.0&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		if PSEUDO_DOUBLE_TOTAL_BITS-exponent>=64  {
			return 0;
		}
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if exponent>PSEUDO_DOUBLE_TOTAL_BITS {
				panic!("Overflow converting PseudoDouble to i64");
			}
		}
		return (x.0&EXP_MASK_INV)>>(PSEUDO_DOUBLE_TOTAL_BITS-exponent);
    }
}

impl From<PseudoDouble> for u64 {
    fn from(x: PseudoDouble) -> Self {
		if x.0<0 {
			panic!("Overflow converting negative PseudoDouble to u64");
		}
		if x.0==0 {
			return 0;
		}
		let exponent=((x.0&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		let vx=(x.0&EXP_MASK_INV) as u64;
		if exponent==PSEUDO_DOUBLE_TOTAL_BITS+1 {
			return vx<<1;
		}
		if PSEUDO_DOUBLE_TOTAL_BITS-exponent>=64  {
			return 0;
		}
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if exponent>PSEUDO_DOUBLE_TOTAL_BITS {
				panic!("Overflow converting PseudoDouble to u64");
			}
		}
		return vx>>(PSEUDO_DOUBLE_TOTAL_BITS-exponent);
    }
}

impl From<PseudoDouble> for f64 {
    fn from(x: PseudoDouble) -> Self {
		if x.0==0 {
			return 0.0;
		}
		let mut vx=x.0&EXP_MASK_INV;
		let sgn;
		let exponent=((x.0&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS+0x3FF-2) as i32;
		if vx<0 {
			sgn=0x8000000000000000u64;
			if vx==(1<<(PSEUDO_DOUBLE_TOTAL_BITS-1)) {
				if exponent< -1 {
					return 0.0;
				}
				if exponent>=0x7FF {
					return f64::NAN;
				}
				return f64::from_bits((((exponent+1) as u64)<<52)+sgn);
			}
			vx=-vx;
		} else {
			sgn=0;
		}
		if exponent<0 {
			return 0.0;
		}
		if exponent>0x7FF {
			return f64::NAN;
		}
		return f64::from_bits(((vx&0x3FFFFFFFFFFFFFFF)>>10) as u64+((exponent as u64)<<52)+sgn);
	}
}

impl From<i32> for PseudoDouble {
    fn from(x : i32) -> Self {
		return PseudoDouble::from(x as i64);
    }
}

impl From<i16> for PseudoDouble {
    fn from(x : i16) -> Self {
		return PseudoDouble::from(x as i64);
    }
}

impl From<i8> for PseudoDouble {
    fn from(x : i8) -> Self {
		return PseudoDouble::from(x as i64);
    }
}

impl From<u32> for PseudoDouble {
    fn from(x : u32) -> Self {
		return PseudoDouble::from(x as u64);
    }
}

impl From<u16> for PseudoDouble {
    fn from(x : u16) -> Self {
		return PseudoDouble::from(x as u64);
    }
}

impl From<u8> for PseudoDouble {
    fn from(x : u8) -> Self {
		return PseudoDouble::from(x as u64);
    }
}

impl From<PseudoDouble> for u32 {
    fn from(x: PseudoDouble) -> Self {
		return u64::from(x) as u32;
	}
}

impl From<PseudoDouble> for u16 {
    fn from(x: PseudoDouble) -> Self {
		return u64::from(x) as u16;
	}
}

impl From<PseudoDouble> for u8 {
    fn from(x: PseudoDouble) -> Self {
		return u64::from(x) as u8;
	}
}

impl From<PseudoDouble> for i32 {
    fn from(x: PseudoDouble) -> Self {
		return i64::from(x) as i32;
	}
}

impl From<PseudoDouble> for i16 {
    fn from(x: PseudoDouble) -> Self {
		return i64::from(x) as i16;
	}
}

impl From<PseudoDouble> for i8 {
    fn from(x: PseudoDouble) -> Self {
		return i64::from(x) as i8;
	}
}

impl Neg for PseudoDouble {
    type Output = Self;
	fn neg(self) -> Self {
		let expx=self.0&EXP_MASK;
		let vx=self.0&EXP_MASK_INV;
		if (vx<<2)==0 {
			let hi_byte=(self.0>>(PSEUDO_DOUBLE_TOTAL_BITS-8)) as u8;
			if hi_byte==0x80 {
				if cfg!(feature="panic_on_pseudodouble_overflow") {
					if expx==EXP_MASK {
						panic!("Overflow in PseudoDouble neg");
					}
				}
				return PseudoDouble(((vx as u64)>>1) as i64+expx+1);
			}
			if hi_byte==0x40 {
				if cfg!(feature="check_on_pseudodouble_underflow") {
					if expx==0 {
						return PD_ZERO;
					}
				}
				return PseudoDouble((vx<<1)+expx-1);
			}
		}
		return PseudoDouble(-vx+expx as i64);
	}
}

impl Add for PseudoDouble {
    type Output = Self;
	fn add(self, other: Self) -> Self {
		let expx=(self.0&EXP_MASK) as i32;
		let expy=(other.0&EXP_MASK) as i32;
		let ydiffx=expy-expx;
		if ydiffx>=(PSEUDO_DOUBLE_TOTAL_BITS-1) {
			return other;
		}
		if ydiffx<=-(PSEUDO_DOUBLE_TOTAL_BITS-1) {
			return self;
		}
		let mut vx=((self.0&EXP_MASK_INV)>>1) as i64;
		let mut vy=((other.0&EXP_MASK_INV)>>1) as i64;
		let exp_max;
		if ydiffx>=0 {
			exp_max=expy+1;
			vx>>=ydiffx;
		} else {
			exp_max=expx+1;
			vy>>=-ydiffx;
		}
		let vr=((vx+vy+PSEUDO_DOUBLE_HALF_ULP)&!PSEUDO_DOUBLE_HALF_ULP) as i64;
		if vr==0 {
			// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
			return PD_ZERO
		} else {
			let mut leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
			if leading_bits>exp_max {
				leading_bits=exp_max;
			}
			let new_exponent=exp_max-leading_bits;
			if cfg!(feature="check_on_pseudodouble_underflow") {
				if new_exponent<0 {
					return PD_ZERO;
				}
			}
			if cfg!(feature="panic_on_pseudodouble_overflow") {
				if new_exponent as u32>EXP_MASK as u32 {
					panic!("Overflow in PseudoDouble add");
				}
			}
			return PseudoDouble(((vr<<leading_bits)&EXP_MASK_INV)+new_exponent as i64);
		}
	}
}

impl Sub for PseudoDouble {
    type Output = Self;
	fn sub(self, other: Self) -> Self {
		let expx=(self.0&EXP_MASK) as i32;
		let expy=(other.0&EXP_MASK) as i32;
		let ydiffx=expy-expx;
		if ydiffx>=(PSEUDO_DOUBLE_TOTAL_BITS-1) {
			return -other;
		}
		if ydiffx<=-(PSEUDO_DOUBLE_TOTAL_BITS-1) {
			return self;
		}
		let mut vx=((self.0&EXP_MASK_INV)>>1) as i64;
		let mut vy=((other.0&EXP_MASK_INV)>>1) as i64;
		let exp_max;
		if ydiffx>=0 {
			exp_max=expy+1;
			vx>>=ydiffx;
		} else {
			exp_max=expx+1;
			vy>>=-ydiffx;
		}
		let vr=((vx-vy+PSEUDO_DOUBLE_HALF_ULP)&!PSEUDO_DOUBLE_HALF_ULP) as i64;
		if vr==0 {
			// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
			return PD_ZERO
		} else {
			let mut leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
			if leading_bits>exp_max {
				leading_bits=exp_max;
			}
			let new_exponent=exp_max-leading_bits;
			if cfg!(feature="check_on_pseudodouble_underflow") {
				if new_exponent<0 {
					return PD_ZERO;
				}
			}
			if cfg!(feature="panic_on_pseudodouble_overflow") {
				if new_exponent as u32>EXP_MASK as u32 {
					panic!("Overflow in PseudoDouble add");
				}
			}
			return PseudoDouble(((vr<<leading_bits)&EXP_MASK_INV)+new_exponent as i64);
		}
	}
}

impl Mul for PseudoDouble {
    type Output = Self;
	fn mul(self, other: Self) -> Self {
		let expx=(self.0&EXP_MASK) as i32;
		let expy=(other.0&EXP_MASK) as i32;
		let vx=(self.0&EXP_MASK_INV) as i128;
		let vy=(other.0&EXP_MASK_INV) as i128;
		let vr=((vx*vy)>>64) as i64;
		if vr==0 {
			return PD_ZERO;
		}
		let leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
		let new_exponent=expx+expy-PSEUDO_DOUBLE_EXP_BIAS as i32-leading_bits;
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if new_exponent as u32>EXP_MASK as u32 {
				panic!("Overflow in PseudoDouble add");
			}
		}
		if cfg!(feature="check_on_pseudodouble_underflow") {
			if new_exponent<0 {
				return PD_ZERO;
			}
		}
		return PseudoDouble(((vr<<leading_bits)&EXP_MASK_INV)+new_exponent as i64);
	}
}

impl Div for PseudoDouble {
    type Output = Self;
	fn div(self, other: Self) -> Self {
		let expx=(self.0&EXP_MASK) as i32;
		let expy=(other.0&EXP_MASK) as i32;
		let vx=(self.0&EXP_MASK_INV) as i128;
		let vy=(other.0&EXP_MASK_INV) as i128;
		if vy==0 { // leave this one in to avoid division by zero signal
			panic!("Division by zero");
		}
		let vxb=((vx>>2) as i128)<<64;
		let vyb=vy as i128;
		let vrb=vxb/vyb;
		let vr=vrb as i64;
		if vr==0 {
			return PD_ZERO;
		}
		let leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
		let new_exponent=2+expx-expy+PSEUDO_DOUBLE_EXP_BIAS as i32-leading_bits;
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if new_exponent as u32>EXP_MASK as u32 {
				panic!("Overflow in PseudoDouble add");
			}
		}
		if cfg!(feature="check_on_pseudodouble_underflow") {
			if new_exponent<0 {
				return PD_ZERO;
			}
		}
		return PseudoDouble(((vr<<leading_bits)&EXP_MASK_INV)+new_exponent as i64);
	}
}

impl AddAssign for PseudoDouble {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

impl SubAssign for PseudoDouble {
    fn sub_assign(&mut self, other: Self) {
        *self = *self - other;
    }
}

impl MulAssign for PseudoDouble {
    fn mul_assign(&mut self, other: Self) {
        *self = *self * other;
    }
}

impl DivAssign for PseudoDouble {
    fn div_assign(&mut self, other: Self) {
        *self = *self / other;
    }
}

impl Ord for PseudoDouble {
    fn cmp(&self, other: &Self) -> Ordering {
		let neg=(self.0>>(PSEUDO_DOUBLE_TOTAL_BITS-1))!=0;
		if ((self.0^other.0)>>(PSEUDO_DOUBLE_TOTAL_BITS-1))!=0 {
			return if neg {Ordering::Less} else {Ordering::Greater};
		}
		// signs are the same, check exponent
		let expdiff=(other.0&EXP_MASK)-(self.0&EXP_MASK);
		if expdiff!=0 {
			return  if (expdiff>0)^neg {Ordering::Less} else {Ordering::Greater};
		} else {
			// exponents are the same so don't need to mask off, check mantissa
			return self.0.cmp(&other.0);
		}
    }
}

impl PartialOrd for PseudoDouble {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
		let neg=(self.0>>(PSEUDO_DOUBLE_TOTAL_BITS-1))!=0;
		if ((self.0^other.0)>>(PSEUDO_DOUBLE_TOTAL_BITS-1))!=0 {
			return if neg {Some(Ordering::Less)} else {Some(Ordering::Greater)};
		}
		// signs are the same, check exponent
		let expdiff=(other.0&EXP_MASK)-(self.0&EXP_MASK);
		if expdiff!=0 {
			return  if (expdiff>0)^neg {Some(Ordering::Less)} else {Some(Ordering::Greater)};
		} else {
			// exponents are the same so don't need to mask off, check mantissa
			return Some(self.0.cmp(&other.0));
		}
    }
}

impl PseudoDouble {

	pub const fn pd_create_mantissa_exp10(dd: i64, ee: i32) -> PseudoDouble {
		if dd==0 {
			return PD_ZERO;
		}
		let mut d=dd;
		let mut e=ee;
		let negative=d<0;
		let mut nexp=0i32;
		while e>0 {
			let lead_bits=(if negative {!d} else {d}).leading_zeros() as i32;
			if lead_bits<5 {
				// check that there is no overflow
				d>>=5-lead_bits;
				nexp+=5-lead_bits;
			}
			d*=10;
			e-=1;
		}
		while e<0 {
			let lead_bits=(if negative {!d} else {d}).leading_zeros() as i32;
			if lead_bits>1 {
				// make the number as accurate as possible
				d<<=lead_bits-1;
				nexp-=lead_bits-1;
			}
			d/=10;
			e+=1;
		}
		let lead_bits=(if negative {!d} else {d}).leading_zeros() as i32;
		let exp=nexp as i64+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits as i64;
		return PseudoDouble(((shift_left_signed(d,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))&EXP_MASK_INV)+exp as i64);
	}

	pub const fn pd_create_mantissa_exp2(d: i64, e: i32) -> PseudoDouble {
		if d==0 {
			return PD_ZERO;
		}
		let lead_bits=(if d<0 {!d} else {d}).leading_zeros() as i32;
		return PseudoDouble(((shift_left_signed(d,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits as i64+e as i64);
	}

	pub const fn gt_zero(x : PseudoDouble) -> bool {
		return x.0>0;
	}

	pub const fn gte_zero(x : PseudoDouble) -> bool {
		return x.0>=0;
	}

	pub const fn lt_zero(x : PseudoDouble) -> bool {
		return x.0<0;
	}

	pub const fn lte_zero(x : PseudoDouble) -> bool {
		return x.0<=0;
	}

	pub fn double_to_pseudodouble(f:f64) -> PseudoDouble {
		if f==0.0 {
			return PD_ZERO;
		}
		let i=f64::to_bits(f) as i64;
		let negative=i<0;
		let raw_exponent=(((i as u64)>>52)&0x7FF) as i64;
		let exponent=raw_exponent+PSEUDO_DOUBLE_EXP_BIAS as i64-0x3FF+2;
		let old_mantissa=i&0xFFFFFFFFFFFFFi64;
		let mantissa=old_mantissa+0x10000000000000i64; // add in the implied bit
		if negative {
			if old_mantissa==0 {
				if exponent<1 {
					return PD_ZERO;
				}
				if exponent>EXP_MASK+1 {
					panic!("Overflow in double_to_pseudodouble");
				}
				return PseudoDouble((1<<(PSEUDO_DOUBLE_TOTAL_BITS-1))+exponent-1);
			}
		}
		if exponent<0 {
			return PD_ZERO;
		}
		if exponent>EXP_MASK {
			panic!("Overflow in double_to_pseudodouble");
		}
		let mantissa=shift_left_signed(mantissa,PSEUDO_DOUBLE_TOTAL_BITS-54);
		//mantissa=(mantissa+PSEUDO_DOUBLE_HALF_ULP)&~PSEUDO_DOUBLE_HALF_ULP;
		if negative {
			return PseudoDouble(-(mantissa&EXP_MASK_INV)+exponent);
		} else {
			return PseudoDouble((mantissa&EXP_MASK_INV)+exponent);
		}
	}

	pub const fn const_neg(self) -> Self {
		let expx=self.0&EXP_MASK;
		let vx=self.0&EXP_MASK_INV;
		if (vx<<2)==0 {
			let hi_byte=(self.0>>(PSEUDO_DOUBLE_TOTAL_BITS-8)) as u8;
			if hi_byte==0x80 {
				if cfg!(feature="panic_on_pseudodouble_overflow") {
					if expx==EXP_MASK {
						panic!("Overflow in PseudoDouble neg");
					}
				}
				return PseudoDouble(((vx as u64)>>1) as i64+expx+1);
			}
			if hi_byte==0x40 {
				if cfg!(feature="check_on_pseudodouble_underflow") {
					if expx==0 {
						return PD_ZERO;
					}
				}
				return PseudoDouble((vx>>1)+expx-1);
			}
		}
		return PseudoDouble(-vx+expx as i64);
	}

	pub const fn const_add(self, other: Self) -> Self {
		let expx=(self.0&EXP_MASK) as i32;
		let expy=(other.0&EXP_MASK) as i32;
		let ydiffx=expy-expx;
		if ydiffx>=(PSEUDO_DOUBLE_TOTAL_BITS-1) {
			return other;
		}
		if ydiffx<=-(PSEUDO_DOUBLE_TOTAL_BITS-1) {
			return self;
		}
		let mut vx=((self.0&EXP_MASK_INV)>>1) as i64;
		let mut vy=((other.0&EXP_MASK_INV)>>1) as i64;
		let exp_max;
		if ydiffx>=0 {
			exp_max=expy+1;
			vx>>=ydiffx;
		} else {
			exp_max=expx+1;
			vy>>=-ydiffx;
		}
		let vr=((vx+vy+PSEUDO_DOUBLE_HALF_ULP)&!PSEUDO_DOUBLE_HALF_ULP) as i64;
		if vr==0 {
			// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
			return PD_ZERO
		} else {
			let mut leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
			if leading_bits>exp_max {
				leading_bits=exp_max;
			}
			let new_exponent=exp_max-leading_bits;
			if cfg!(feature="check_on_pseudodouble_underflow") {
				if new_exponent<0 {
					return PD_ZERO;
				}
			}
			if cfg!(feature="panic_on_pseudodouble_overflow") {
				if new_exponent as u32>EXP_MASK as u32 {
					panic!("Overflow in PseudoDouble add");
				}
			}
			return PseudoDouble(((vr<<leading_bits)&EXP_MASK_INV)+new_exponent as i64);
		}
	}

	pub const fn const_sub(self, other: Self) -> Self {
		let expx=(self.0&EXP_MASK) as i32;
		let expy=(other.0&EXP_MASK) as i32;
		let ydiffx=expy-expx;
		if ydiffx>=(PSEUDO_DOUBLE_TOTAL_BITS-1) {
			return other.const_neg();
		}
		if ydiffx<=-(PSEUDO_DOUBLE_TOTAL_BITS-1) {
			return self;
		}
		let mut vx=((self.0&EXP_MASK_INV)>>1) as i64;
		let mut vy=((other.0&EXP_MASK_INV)>>1) as i64;
		let exp_max;
		if ydiffx>=0 {
			exp_max=expy+1;
			vx>>=ydiffx;
		} else {
			exp_max=expx+1;
			vy>>=-ydiffx;
		}
		let vr=((vx-vy+PSEUDO_DOUBLE_HALF_ULP)&!PSEUDO_DOUBLE_HALF_ULP) as i64;
		if vr==0 {
			// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
			return PD_ZERO
		} else {
			let mut leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
			if leading_bits>exp_max {
				leading_bits=exp_max;
			}
			let new_exponent=exp_max-leading_bits;
			if cfg!(feature="check_on_pseudodouble_underflow") {
				if new_exponent<0 {
					return PD_ZERO;
				}
			}
			if cfg!(feature="panic_on_pseudodouble_overflow") {
				if new_exponent as u32>EXP_MASK as u32 {
					panic!("Overflow in PseudoDouble add");
				}
			}
			return PseudoDouble(((vr<<leading_bits)&EXP_MASK_INV)+new_exponent as i64);
		}
	}

	pub const fn const_mul(self, other: Self) -> Self {
		let expx=(self.0&EXP_MASK) as i32;
		let expy=(other.0&EXP_MASK) as i32;
		let vx=(self.0&EXP_MASK_INV) as i128;
		let vy=(other.0&EXP_MASK_INV) as i128;
		let vr=((vx*vy)>>64) as i64;
		if vr==0 {
			return PD_ZERO;
		}
		let leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
		let new_exponent=expx+expy-PSEUDO_DOUBLE_EXP_BIAS as i32-leading_bits;
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if new_exponent as u32>EXP_MASK as u32 {
				panic!("Overflow in PseudoDouble add");
			}
		}
		if cfg!(feature="check_on_pseudodouble_underflow") {
			if new_exponent<0 {
				return PD_ZERO;
			}
		}
		return PseudoDouble(((vr<<leading_bits)&EXP_MASK_INV)+new_exponent as i64);
	}

	pub const fn const_div(self, other: Self) -> Self {
		let expx=(self.0&EXP_MASK) as i32;
		let expy=(other.0&EXP_MASK) as i32;
		let vx=(self.0&EXP_MASK_INV) as i128;
		let vy=(other.0&EXP_MASK_INV) as i128;
		if vy==0 { // leave this one in to avoid division by zero signal
			panic!("Division by zero");
		}
		let vr=((vx>>2) as i128/vy) as i64;
		if vr==0 {
			return PD_ZERO;
		}
		let leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
		let new_exponent=2+expx-expy+PSEUDO_DOUBLE_EXP_BIAS as i32-leading_bits;
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if new_exponent as u32>EXP_MASK as u32 {
				panic!("Overflow in PseudoDouble add");
			}
		}
		if cfg!(feature="check_on_pseudodouble_underflow") {
			if new_exponent<0 {
				return PD_ZERO;
			}
		}
		return PseudoDouble((vr&EXP_MASK_INV)+new_exponent as i64);
	}

    pub const fn const_less_than_or_equal(&self, other: &Self) -> bool {
		if self.0==other.0 {
			return true;
		}
		let neg=(self.0>>(PSEUDO_DOUBLE_TOTAL_BITS-1))!=0;
		if ((self.0^other.0)>>(PSEUDO_DOUBLE_TOTAL_BITS-1))!=0 {
			return neg;
		}
		// signs are the same, check exponent
		let expdiff=(other.0&EXP_MASK)-(self.0&EXP_MASK);
		if expdiff!=0 {
			return (expdiff>0)^neg;
		} else {
			// exponents are the same so don't need to mask off, check mantissa
			return self.0<other.0;
		}
    }

    pub const fn const_less_than(&self, other: &Self) -> bool {
		if self.0==other.0 {
			return false;
		}
		let neg=(self.0>>(PSEUDO_DOUBLE_TOTAL_BITS-1))!=0;
		if ((self.0^other.0)>>(PSEUDO_DOUBLE_TOTAL_BITS-1))==0 {
			return neg;
		}
		// signs are the same, check exponent
		let expdiff=(other.0&EXP_MASK)-(self.0&EXP_MASK);
		if expdiff!=0 {
			return (expdiff>0)^neg;
		} else {
			// exponents are the same so don't need to mask off, check mantissa
			return self.0<other.0;
		}
    }

    pub const fn floor(self) -> PseudoDouble {
		let exponent=self.0&EXP_MASK;
		let e=(exponent-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		if e<2 {
			return if self.0<0 {PD_NEG_ONE} else {PD_ZERO};
		}
		if e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS {
			return self;
		}
		let m=(1<<(PSEUDO_DOUBLE_TOTAL_BITS-e))-1;
		return PseudoDouble((self.0&!m)+exponent);
    }

	pub const fn ceil(self) -> PseudoDouble {
		let exponent=self.0&EXP_MASK;
		let e=(exponent-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		let vx=self.0&EXP_MASK_INV;
		if e<2 {
			if e==1 && (vx<<1)==0 { // special test for ceil(-1)=-1
				return self;
			}
			return if self.0>0 {PD_ONE} else {PD_ZERO};
		}
		if e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS {
			return self;
		}
		let m=(1<<(PSEUDO_DOUBLE_TOTAL_BITS-e-1))-1;
		let vr=((vx>>1)+m)&!m;
		let leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i64 - 1;
		let new_exponent=exponent+1-leading_bits;
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if new_exponent as u32>EXP_MASK as u32 {
				panic!("Overflow in PseudoDouble add");
			}
		}
		if cfg!(feature="check_on_pseudodouble_underflow") {
			if new_exponent<0 {
				return PD_ZERO;
			}
		}
		return PseudoDouble((vr<<leading_bits)+new_exponent);
    }

	pub const fn round(self) -> PseudoDouble {
		let exponent=self.0&EXP_MASK;
		let e=(exponent-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		let vx=self.0&EXP_MASK_INV;
		if e<1 {
			if e==0 && (vx<<1)==0 { // special test for round(-0.5)=-1
				return PseudoDouble(self.0+1);
			}
			return PD_ZERO;
		}
		if e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS {
			return self;
		}
		let add=1<<(PSEUDO_DOUBLE_TOTAL_BITS-e-2);
		let m=(add<<1)-1;
		let vr=((vx>>1)+if vx>0 {add} else {add-1})&!m;
		let leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i64 - 1;
		let new_exponent=exponent+1-leading_bits;
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if new_exponent as u32>EXP_MASK as u32 {
				panic!("Overflow in PseudoDouble add");
			}
		}
		if cfg!(feature="check_on_pseudodouble_underflow") {
			if new_exponent<0 {
				return PD_ZERO;
			}
		}
		return PseudoDouble((vr<<leading_bits)+new_exponent);
    }

	pub const fn abs(self) -> PseudoDouble {
		if self.0>=0 {
			return self;
		}
		let expx=self.0&EXP_MASK;
		let vx=self.0&EXP_MASK_INV;
		if (vx<<2)==0 {
			let hi_byte=(self.0>>(PSEUDO_DOUBLE_TOTAL_BITS-8)) as u8;
			if hi_byte==0x80 {
				if cfg!(feature="panic_on_pseudodouble_overflow") {
					if expx==EXP_MASK {
						panic!("Overflow in PseudoDouble abs");
					}
				}
				return PseudoDouble((vx>>1)+expx+1);
			}
		}
		return PseudoDouble(-vx+expx as i64);
	}

	pub const fn inv_sqrt(self) -> PseudoDouble {
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if self.0<0 {
				panic!("sqrt of negative number");
			}
		}
		let mut exponent=self.0&EXP_MASK;
		let mut mantissa=self.0&EXP_MASK_INV;
		// [01.00 .. 11.11] = [2^0 .. 2^2)  ->  [2^0 .. 2^-1) = [1 .. 0.5)
		if (exponent&1)!=0 {
			exponent-=1;
			mantissa<<=1;
		} else {
			if (mantissa<<2)==0 {
				return PseudoDouble(mantissa+3*(PSEUDO_DOUBLE_EXP_BIAS>>1)+3-(exponent>>1));
			}
		}
		return PseudoDouble((inv_sqrt64_fixed(mantissa as u64) as i64&EXP_MASK_INV)+3*(PSEUDO_DOUBLE_EXP_BIAS>>1)+2-(exponent>>1));
	}

	pub const fn sqrt(self) -> PseudoDouble {
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if self.0<0 {
				panic!("sqrt of negative number");
			}
		}
		if self.0==0 {
			return self;
		}
		let mut exponent=self.0&EXP_MASK;
		let mut mantissa=self.0&EXP_MASK_INV;
		if (exponent&1)!=0 {
			exponent-=1;
			mantissa<<=1;
		} else {
			if (mantissa<<2)==0 {
				return PseudoDouble(mantissa+(PSEUDO_DOUBLE_EXP_BIAS>>1)+1+(exponent>>1));
			}
		}
		// (1,4) * (1,0.5) = (1,2)
		let y=(multu64hi(inv_sqrt64_fixed(mantissa as u64>>(64-PSEUDO_DOUBLE_TOTAL_BITS))<<(64-PSEUDO_DOUBLE_TOTAL_BITS),mantissa as u64)<<1) as i64;
		return PseudoDouble((y&EXP_MASK_INV)+(PSEUDO_DOUBLE_EXP_BIAS>>1)+1+(exponent>>1));
	}

	pub const fn ldexp(self, y:i32) -> PseudoDouble {
		let yy=y as i64;
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if (self.0&EXP_MASK)+yy>EXP_MASK {
				panic!("Overflow in PseudoDouble ldexp");
			}
		}
		if cfg!(feature="check_on_pseudodouble_underflow") {
			if (self.0&EXP_MASK)+yy<0 {
				return PD_ZERO;
			}
		}
		return PseudoDouble(self.0+yy);
	}

	// e^x=(2^log2(e))^x=2^(log2(e)*x)
	pub const fn exp2(self) -> PseudoDouble {
		if self.0==0 {
			return PD_ONE;
		}
		let exponent=self.0&EXP_MASK;
		let e=(exponent-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		let new_exponent;
		let mut fraction;
		if e<2 {
			if self.0<0 {
				new_exponent=-1;
				if e==1 {
					fraction=(self.0&EXP_MASK_INV)<<1;
				} else {
					fraction=(self.0&EXP_MASK_INV)>>-e;
				}
			} else {
				new_exponent=0;
				fraction=self.0&EXP_MASK_INV;
				if e==1 {
					fraction<<=1;
				} else {
					fraction>>=-e;
				}
			}
		} else if e<=PSEUDO_DOUBLE_EXP_BITS { // maself.0=2^(2^PSEUDO_DOUBLE_EXP_BITS)), log2(maself.0)=2^PSEUDO_DOUBLE_EXP_BITS
			let m=(1i64<<(PSEUDO_DOUBLE_TOTAL_BITS-e))-1;
			new_exponent=(self.0&!m)>>(PSEUDO_DOUBLE_TOTAL_BITS-e);
			fraction=(self.0&EXP_MASK_INV&m)<<e;
		} else {
			// common to have underflow, so leave this in even if errors otherwise are turned off
			if self.0<0 {
				return PD_ZERO;
			} else {
				if cfg!(feature="panic_on_pseudodouble_overflow") {
					panic!("Overflow in PseudoDouble exp");
				} else {
					// invalid result, but keep the compiler happy
					new_exponent=0;
					fraction=0;
				}
			}
		}
		let newe=new_exponent+PSEUDO_DOUBLE_EXP_BIAS+2;
		if newe<0 {
			return PD_ZERO;
		} else {
			if cfg!(feature="panic_on_pseudodouble_overflow") {
				if newe>EXP_MASK {
					panic!("Overflow in PseudoDouble exp");
				}
			}
		}
		return PseudoDouble(newe+(((exp2_64_fixed((fraction<<(64-PSEUDO_DOUBLE_TOTAL_BITS)) as u64) as i64)<<(64-PSEUDO_DOUBLE_TOTAL_BITS))&EXP_MASK_INV));
	}

	pub const fn log2(self) -> PseudoDouble {
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if self.0<=0 {
				panic!("PseudoDouble log2 of non-positive number");
			}
		}
		let exponent=self.0&EXP_MASK;
		let e=exponent-PSEUDO_DOUBLE_EXP_BIAS-2;
		let mantissa=((self.0&EXP_MASK_INV)<<2) as u64>>1;
		let log_frac=log2_64_fixed(mantissa);
		if e==0 {
			if log_frac==0 {
				return PD_ZERO;
			}
			let lead_bits=log_frac.leading_zeros() as i64;
			return PseudoDouble((((log_frac as i64)<<(lead_bits-1))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+2-lead_bits);
		} else if e==-1 {
			let logfrac2=log_frac+0x8000000000000000u64;
			let lead_bits=(!logfrac2).leading_zeros() as i64;
			return PseudoDouble((((logfrac2 as i64)<<(lead_bits-1))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+2-lead_bits);
		}
		let lead_bits=(if e<0 {!e} else {e}).leading_zeros() as i64;
		return PseudoDouble((((e<<(PSEUDO_DOUBLE_TOTAL_BITS+(lead_bits as i32)-65))+((log_frac as i64)>>(64-(lead_bits as i32))))&EXP_MASK_INV) as i64+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits);
	}

// pseudo_double_i pdi_log2(pseudo_double_i x) {
// #if PD_ERROR_CHECK
// 	if(((signed_pd_internal)x)<=0) {
// 		PD_DO_ERROR_RANGE;
// 	}
// #endif
// 	int64_t exponent=(x&EXP_MASK);
// 	int64_t e=exponent-PSEUDO_DOUBLE_EXP_BIAS-2;
// 	uint64_t mantissa=((x&EXP_MASK_INV)<<2)>>1;
// 	uint64_t log_frac=log2_64_fixed(mantissa);
// 	if(e==0) {
// 		if(log_frac==0) {
// 			return 0;
// 		}
// 		int lead_bits=clz(log_frac);
// 		return ((log_frac<<(lead_bits-1))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+2-lead_bits;
// 	} else if(e==-1) {
// 		log_frac+=0x8000000000000000ULL;
// 		int lead_bits=clz(~log_frac);
// 		return ((log_frac<<(lead_bits-1))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+2-lead_bits;
// 	}
// 	int negative=(e<0);
// 	int lead_bits=clz(negative?~e:e);
// 	return (((e<<(PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))+(log_frac>>(64-lead_bits)))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits;
// }

	// x^y = e^ln(x)^y = 2^(y*ln2(x))
	pub const fn pow(self, y:Self) -> Self {
		if cfg!(feature="panic_on_pseudodouble_overflow") {
			if self.0<=0 {
				panic!("PseudoDouble pow of non-positive number");
			}
		}
		let exponent=self.0&EXP_MASK;
		let e=exponent-PSEUDO_DOUBLE_EXP_BIAS-2;
		let mantissa=((self.0&EXP_MASK_INV)<<2) as u64>>1;
		let mut log_frac=log2_64_fixed(mantissa);
		let vx;
		let expx;
		if e==0 {
			if log_frac==0 {
				return PD_ONE;
			}
			let lead_bits=log_frac.leading_zeros() as i32;
			vx=(log_frac<<(lead_bits-1)) as i64;
			expx=2-lead_bits;
		} else if e==-1 {
			log_frac+=0x8000000000000000u64;
			let lead_bits=(!log_frac).leading_zeros() as i32;
			vx=(log_frac<<(lead_bits-1)) as i64;
			expx=2-lead_bits;
		} else {
			let lead_bits=(if e<0 {!e} else {e}).leading_zeros() as i32;
			vx=(e<<(PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))+((log_frac as i64)>>(64-lead_bits));
			expx=65-lead_bits;
		}
		let expy=((y.0&EXP_MASK)-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		let vy=y.0&EXP_MASK_INV;
		let mut vr=mults64hi(vx,vy);
		if vr==0 {
			// special case - a mantissa of zero will always make the whole word zero. Makes comparisons much easier
			return PseudoDouble((1i64<<(PSEUDO_DOUBLE_TOTAL_BITS-2))+PSEUDO_DOUBLE_EXP_BIAS+2); // 2^0=1
		}
		let leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
		vr<<=leading_bits;
		let er=expx+expy-leading_bits;
		let new_exponent;
		let mut fraction;
		if er<2 {
			if vr<0 {
				new_exponent=-1;
				if er==1 {
					fraction=vr<<1;
				} else {
					fraction=vr>>-er;
				}
			} else {
				new_exponent=0;
				fraction=vr;
				if er==1 {
					fraction<<=1;
				} else {
					fraction>>=-er;
				}
			}
		} else if er<=PSEUDO_DOUBLE_EXP_BITS { // max=2^(2^PSEUDO_DOUBLE_EXP_BITS)), log2(max)=2^PSEUDO_DOUBLE_EXP_BITS
			let m=(1i64<<(PSEUDO_DOUBLE_TOTAL_BITS-er))-1;
			new_exponent=((vr&!m))>>(PSEUDO_DOUBLE_TOTAL_BITS-er);
			fraction=(vr&m)<<er;
		} else {
			if vr<0 {
				return PD_ZERO;
			} else {
				if cfg!(feature="panic_on_pseudodouble_overflow") {
					panic!("PseudoDouble overflow on pow");
				} else {
					// invalid result, but keep the compiler happy
					new_exponent=0;
					fraction=0;
				}
			}
		}
		let newe=new_exponent+PSEUDO_DOUBLE_EXP_BIAS+2;
		if newe<0 { // common to have underflow, so leave this in
			return PD_ZERO;
		} else {
			if cfg!(feature="panic_on_pseudodouble_overflow") {
				if newe>EXP_MASK {
					panic!("PseudoDouble overflow on pow");
				}
			}
		}
		return PseudoDouble(newe+(((exp2_64_fixed((fraction<<(64-PSEUDO_DOUBLE_TOTAL_BITS)) as u64) as i64)<<(64-PSEUDO_DOUBLE_TOTAL_BITS))&EXP_MASK_INV));
	}

	pub const fn exp(self) -> PseudoDouble {
		return self.const_mul(PD_LOG_2_E).exp2();
	}

	pub const fn exp10(self) -> PseudoDouble {
		return self.const_mul(PD_LOG_2_10).exp2();
	}

	pub const fn log(self) -> PseudoDouble {
		return self.log2().const_mul(PD_INV_LOG_2_E);
	}

	pub const fn log10(self) -> PseudoDouble {
		return self.log2().const_mul(PD_INV_LOG_2_10);
	}

	pub fn string_to_pdi(s:&str) -> Option<PseudoDouble> {
		let mut chars = s.bytes();
		let neg;
		let mut ch=chars.next();
		if ch==None {
			return None;
		}
		if ch.unwrap()==b'-' {
			neg=true;
			ch=chars.next();
			if ch==None {
				return None;
			}
		} else {
			neg=false;
		}
		let mut acc=0i64;
		while ch!=None && ch.unwrap()>=b'0' && (ch.unwrap()<=b'9') {
			acc=acc*10+(ch.unwrap()-b'0') as i64;
			ch=chars.next();
		}
		if ch==None {
			return Some(PseudoDouble::from(if neg {-acc} else {acc}));
		}
		if ch.unwrap()!=b'.' {
			return None;
		}
		ch=chars.next();
		let mut frac_num=0i64;
		let mut frac_den=1i64;
		while ch!=None && ch.unwrap()>=b'0' && (ch.unwrap()<=b'9') {
			frac_num=frac_num*10+(ch.unwrap()-b'0') as i64;
			frac_den*=10;
			ch=chars.next();
		}
		if ch!=None {
			return None;
		}
		let ret=PseudoDouble::from(acc)+PseudoDouble::from(frac_num)/PseudoDouble::from(frac_den);
		if neg {
			return Some(-ret);
		} else {
			return Some(ret);
		}
	}

	pub const fn sin_rev(self) -> PseudoDouble {
		if self.0==0 {
			return PD_ZERO;
		}
		let exponent=self.0&EXP_MASK;
		let e=(exponent-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		let mut fraction;
		if e<2 {
			if e==1 {
				fraction=(self.0&EXP_MASK_INV)<<1;
			} else {
				fraction=(self.0&EXP_MASK_INV)>>-e;
			}
		} else {
			let m=(1i64<<(PSEUDO_DOUBLE_TOTAL_BITS-e))-1;
			fraction=(self.0&EXP_MASK_INV&m)<<e;
		}
		let negative=fraction<0;
		if negative {
			fraction=-fraction;
		}
		let mut ufraction=fraction as u64;
		if (ufraction>>62)!=0 {
			ufraction=0x8000000000000000u64-ufraction;
		}
		let mut d=sin_rev_64_fixed(ufraction as u64) as i64;
		if d==0 {
			return PD_ZERO;
		}
		if negative {
			d=-d;
		}
		let lead_bits=(if negative {!d} else {d}).leading_zeros() as i32;
		return PseudoDouble((shift_left_signed(d,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65)&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+3-lead_bits as i64);
	}

	pub const fn cos_rev(self) -> PseudoDouble {
		if self.0==0 {
			return PD_ONE;
		}
		let exponent=self.0&EXP_MASK;
		let e=(exponent-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		let mut fraction;
		if e<2 {
			if e==1 {
				fraction=(self.0&EXP_MASK_INV)<<1;
			} else {
				fraction=(self.0&EXP_MASK_INV)>>-e;
			}
		} else {
			let m=(1i64<<(PSEUDO_DOUBLE_TOTAL_BITS-e))-1;
			fraction=(self.0&EXP_MASK_INV&m)<<e;
		}
		if fraction<0x4000000000000000i64 { // _only_ line that is different between sin and cos
			fraction+=0x4000000000000000i64;
		} else {
			fraction=((fraction as u64)+0x4000000000000000u64) as i64;
		}
		let negative=fraction<0;
		if negative {
			fraction=-fraction;
		}
		let mut ufraction=fraction as u64;
		if (ufraction>>62)!=0 {
			ufraction=0x8000000000000000u64-ufraction;
		}
		let mut d=sin_rev_64_fixed(ufraction) as i64;
		if d==0 {
			return PD_ZERO;
		}
		if negative {
			d=-d;
		}
		let lead_bits=(if negative {!d} else {d}).leading_zeros() as i32;
		return PseudoDouble((shift_left_signed(d,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65)&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+3-lead_bits as i64);
	}

	pub const fn atan2_rev(self, other: PseudoDouble) -> PseudoDouble {
		let mut negative=false;
		let add_const;
		let mut y=self;
		let mut x=other;
		if y.0==0 {
			if x.0>=0 {
				return PD_ZERO;
			} else {
				return Self::pd_create_mantissa_exp2(1,-1); // 1/2
			}
		} else if y.0>0 {
			if x.0==0 {
				return Self::pd_create_mantissa_exp2(1,-2); // 1/4
			} else if x.0>0 {
				if y.const_less_than_or_equal(&x) {
					// q1
					add_const=0;
				} else {
					// q2
					{ let t=y; y=x; x=t;}
					add_const=0x4000000000000000u64;
					negative=true;
				}
			} else { // x<0
				x=x.const_neg();
				if y.const_less_than_or_equal(&x) {
					// q4
					add_const=0x8000000000000000u64;
					negative=true;
				} else {
					// q3
					{ let t=y; y=x; x=t;}
					add_const=0x4000000000000000u64;
				}
			}
		} else { // y<0
			y=y.const_neg();
			if x.0==0 {
				return Self::pd_create_mantissa_exp2(3,-2); // 3/2
			} else if x.0>0 {
				if y.const_less_than_or_equal(&x) {
					// q8
					add_const=0; // boolean
					negative=true;
				} else {
					// q7
					{ let t=y; y=x; x=t;}
					add_const=0xC000000000000000u64;
				}
			} else { // x<0
				x=x.const_neg();
				if y.const_less_than_or_equal(&x) {
					// q5
					add_const=0x8000000000000000u64;
				} else {
					// q6
					{ let t=y; y=x; x=t;}
					add_const=0xC000000000000000u64;
					negative=true; // boolean
				}
			}
		}
		let expx=(x.0&EXP_MASK) as i32;
		let expy=(y.0&EXP_MASK) as i32;
		let vx=x.0&EXP_MASK_INV;
		let vy=y.0&EXP_MASK_INV;
		let ratio;
		if x.0==y.0 {
			ratio=0x4000000000000000i64;
		} else if x.0==y.const_neg().0 {
			ratio=-0x4000000000000000i64;
		} else {
			let mut vr=divs64hi(vy>>2,vx);
			if vr==0 {
				ratio=0;
			} else {
				let leading_bits=vr.leading_zeros() as i32-1;
				vr<<=leading_bits;
				let new_exponent=expy-expx-leading_bits;
				if new_exponent< -63 {
					ratio=0;
				} else {
					ratio=vr>>(-new_exponent);
				}
			}
		}
		let mut d=(atan_rev_64_fixed(ratio as u64)>>1) as i128;
		d=add_const as i128+if negative {-d} else {d};
		let d64=d as i64;
		if d64==0 {
			return PD_ZERO;
		}
		let negatived=d<0;
		let lead_bits=(if negatived {!d64} else {d64}).leading_zeros() as i32;
		return PseudoDouble(((shift_left_signed(d64,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+1-lead_bits as i64);
	}

	pub const fn sin(self) -> PseudoDouble {
		 self.const_mul(PD_INV_TAU).sin_rev()
	}

	pub const fn cos(self) -> PseudoDouble {
		 self.const_mul(PD_INV_TAU).cos_rev()
	}

	pub const fn atan2(self, other: PseudoDouble) -> PseudoDouble {
		self.atan2_rev(other).const_mul(PD_TAU)
	}
}
