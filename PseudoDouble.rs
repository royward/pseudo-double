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

#[derive(Debug, Copy, Clone, Eq, PartialEq, PartialOrd)]
pub struct PseudoDouble(i64);

const PSEUDO_DOUBLE_TOTAL_BITS: i32 = 64;
const PSEUDO_DOUBLE_EXP_BITS: i32 = 16;
const EXP_MASK: i64 = (1<<PSEUDO_DOUBLE_EXP_BITS)-1;
const EXP_MASK_INV: i64 = !EXP_MASK as i64;
const PSEUDO_DOUBLE_HALF_ULP: i64 = (1<<(PSEUDO_DOUBLE_EXP_BITS-1))-1;
const PSEUDO_DOUBLE_EXP_BIAS: i64 = 1<<(PSEUDO_DOUBLE_EXP_BITS-1);

#[inline]
fn shift_left_signed(x:i64, shift:i32) -> i64 {
	if shift>0 {x<<shift} else {x>>-shift}
}

#[inline]
fn shift_left_unsigned(x:u64, shift:i32) -> u64 {
	if shift>0 {x<<shift} else {x>>-shift}
}

#[inline]
fn multu64hi(x:u64, y:u64) -> u64 {
	((x as u128*y as u128)>>64) as u64
}

#[inline]
fn mults64hi(x:i64, y:i64) -> i64 {
	((x as i128*y as i128)>>64) as i64
}

// x is a 2.62 unsigned fixed in the range (1,4)
// result is 1.63 unsigned fixed in the range (0.5,1)
fn inv_sqrt64_fixed(x:u64) -> u64 {
	// start with a linear interpolation correct at the endpoints
	// 7/6 - 1/6 x, so 1->1, 4->0.5
	let mut y=3074457345618258602u64-multu64hi(x,12297829382473034410u64);
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
fn exp2_64_fixed(x:u64) -> u64 {
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
fn log2_64_fixed(xu:u64) -> u64 {
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

impl From<i64> for PseudoDouble {
    fn from(x : i64) -> Self {
		if x==0 {
			return PseudoDouble(0);
		} else {
			let lead_bits=(if x<0 {!x} else {x}).leading_zeros() as i32;
			return PseudoDouble(((shift_left_signed(x,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits as i64);
		}
    }
}

impl From<u64> for PseudoDouble {
    fn from(x : u64) -> Self {
		if x==0 {
			return PseudoDouble(0);
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
		if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
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
		if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
			if exponent>PSEUDO_DOUBLE_TOTAL_BITS {
				panic!("Overflow converting PseudoDouble to u64");
			}
		}
		return vx>>(PSEUDO_DOUBLE_TOTAL_BITS-exponent);
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
				if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
					if expx==EXP_MASK {
						panic!("Overflow in PseudoDouble neg");
					}
				}
				return PseudoDouble((vx>>1)+expx+1);
			}
			if hi_byte==0x40 {
				if cfg!(CHECK_ON_PSEUDODOUBLE_UNDERFLOW) {
					if expx==0 {
						return PseudoDouble(0);
					}
				}
				return PseudoDouble((vx>>1)+expx+1);
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
			return PseudoDouble(0)
		} else {
			let mut leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
			if leading_bits>exp_max {
				leading_bits=exp_max;
			}
			let new_exponent=exp_max-leading_bits;
			if cfg!(CHECK_ON_PSEUDODOUBLE_UNDERFLOW) {
				if new_exponent<0 {
					return PseudoDouble(0);
				}
			}
			if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
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
			return PseudoDouble(0)
		} else {
			let mut leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
			if leading_bits>exp_max {
				leading_bits=exp_max;
			}
			let new_exponent=exp_max-leading_bits;
			if cfg!(CHECK_ON_PSEUDODOUBLE_UNDERFLOW) {
				if new_exponent<0 {
					return PseudoDouble(0);
				}
			}
			if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
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
			return PseudoDouble(0);
		}
		let leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
		let new_exponent=expx+expy-PSEUDO_DOUBLE_EXP_BIAS as i32-leading_bits;
		if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
			if new_exponent as u32>EXP_MASK as u32 {
				panic!("Overflow in PseudoDouble add");
			}
		}
		if cfg!(CHECK_ON_PSEUDODOUBLE_UNDERFLOW) {
			if new_exponent<0 {
				return PseudoDouble(0);
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
		let vr=((vx>>2) as i128/vy) as i64;
		if vr==0 {
			return PseudoDouble(0);
		}
		let leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i32 - 1;
		let new_exponent=2+expx-expy+PSEUDO_DOUBLE_EXP_BIAS as i32-leading_bits;
		if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
			if new_exponent as u32>EXP_MASK as u32 {
				panic!("Overflow in PseudoDouble add");
			}
		}
		if cfg!(CHECK_ON_PSEUDODOUBLE_UNDERFLOW) {
			if new_exponent<0 {
				return PseudoDouble(0);
			}
		}
		return PseudoDouble((vr&EXP_MASK_INV)+new_exponent as i64);
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
		if ((self.0^other.0)>>(PSEUDO_DOUBLE_TOTAL_BITS-1))==0 {
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

pub fn gt_zero(x : PseudoDouble) -> bool {
	return x.0>0;
}

pub fn gte_zero(x : PseudoDouble) -> bool {
	return x.0>=0;
}

pub fn lt_zero(x : PseudoDouble) -> bool {
	return x.0<0;
}

pub fn lte_zero(x : PseudoDouble) -> bool {
	return x.0<=0;
}

impl PseudoDouble {

	pub fn floor(self) -> PseudoDouble {
		let exponent=self.0&EXP_MASK;
		let e=(exponent-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		if e<2 {
			return if self.0<0 {PseudoDouble(-1)} else {PseudoDouble(0)};
		}
		if e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS {
			return self;
		}
		let m=(1<<(PSEUDO_DOUBLE_TOTAL_BITS-e))-1;
		return PseudoDouble((self.0&!m)+exponent);
    }

	pub fn ceil(self) -> PseudoDouble {
		let exponent=self.0&EXP_MASK;
		let e=(exponent-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		let vx=self.0&EXP_MASK_INV;
		if e<2 {
			if e==1 && (vx<<1)==0 { // special test for ceil(-1)=-1
				return self;
			}
			return if self.0>0 {PseudoDouble(1)} else {PseudoDouble(0)};
		}
		if e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS {
			return self;
		}
		let m=(1<<(PSEUDO_DOUBLE_TOTAL_BITS-e-1))-1;
		let mut vr=((vx>>1)+m)&!m;
		let leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i64 - 1;
		vr<<=leading_bits;
		let new_exponent=exponent+1-leading_bits;
		if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
			if new_exponent as u32>EXP_MASK as u32 {
				panic!("Overflow in PseudoDouble add");
			}
		}
		if cfg!(CHECK_ON_PSEUDODOUBLE_UNDERFLOW) {
			if new_exponent<0 {
				return PseudoDouble(0);
			}
		}
		return PseudoDouble((vr<<leading_bits)+new_exponent);
    }

	pub fn round(self) -> PseudoDouble {
		let exponent=self.0&EXP_MASK;
		let e=(exponent-PSEUDO_DOUBLE_EXP_BIAS) as i32;
		let vx=self.0&EXP_MASK_INV;
		if e<1 {
			if e==0 && (vx<<1)==0 { // special test for round(-0.5)=-1
				return PseudoDouble(self.0+1);
			}
			return PseudoDouble(0);
		}
		if e>=PSEUDO_DOUBLE_TOTAL_BITS-PSEUDO_DOUBLE_EXP_BITS {
			return self;
		}
		let add=1<<(PSEUDO_DOUBLE_TOTAL_BITS-e-2);
		let m=(add<<1)-1;
		let mut vr=((vx>>1)+if vx>0 {add} else {add-1})&!m;
		let leading_bits=(if vr>0 {vr} else {!vr}).leading_zeros() as i64 - 1;
		vr<<=leading_bits;
		let new_exponent=exponent+1-leading_bits;
		if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
			if new_exponent as u32>EXP_MASK as u32 {
				panic!("Overflow in PseudoDouble add");
			}
		}
		if cfg!(CHECK_ON_PSEUDODOUBLE_UNDERFLOW) {
			if new_exponent<0 {
				return PseudoDouble(0);
			}
		}
		return PseudoDouble((vr<<leading_bits)+new_exponent);
    }

	pub fn abs(self) -> PseudoDouble {
		if self.0>=0 {
			return self;
		}
		let expx=self.0&EXP_MASK;
		let vx=self.0&EXP_MASK_INV;
		if (vx<<2)==0 {
			let hi_byte=(self.0>>(PSEUDO_DOUBLE_TOTAL_BITS-8)) as u8;
			if hi_byte==0x80 {
				if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
					if expx==EXP_MASK {
						panic!("Overflow in PseudoDouble abs");
					}
				}
				return PseudoDouble((vx>>1)+expx+1);
			}
		}
		return PseudoDouble(-vx+expx as i64);
	}

	pub fn inv_sqrt(self) -> PseudoDouble {
		if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
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

	pub fn sqrt(self) -> PseudoDouble {
		if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
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

	pub fn ldexp(self, y:i32) -> PseudoDouble {
		let yy=y as i64;
		if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
			if (self.0&EXP_MASK)+yy>EXP_MASK {
				panic!("Overflow in PseudoDouble ldexp");
			}
		}
		if cfg!(CHECK_ON_PSEUDODOUBLE_UNDERFLOW) {
			if (self.0&EXP_MASK)+yy<0 {
				return PseudoDouble(0);
			}
		}
		return PseudoDouble(self.0+yy);
	}

	// e^x=(2^log2(e))^x=2^(log2(e)*x)
	pub fn exp2(self) -> PseudoDouble {
		if self.0==0 {
			return PseudoDouble((1i64<<(PSEUDO_DOUBLE_TOTAL_BITS-2))+PSEUDO_DOUBLE_EXP_BIAS+2);
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
				return PseudoDouble(0);
			} else {
				if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
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
			return PseudoDouble(0);
		} else {
			if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
				if newe>EXP_MASK {
					panic!("Overflow in PseudoDouble exp");
				}
			}
		}
		return PseudoDouble(newe+(((exp2_64_fixed((fraction<<(64-PSEUDO_DOUBLE_TOTAL_BITS)) as u64) as i64)<<(64-PSEUDO_DOUBLE_TOTAL_BITS))&EXP_MASK_INV));
	}

	pub fn log2(self) -> PseudoDouble {
		if cfg!(PANIC_ON_PSEUDODOUBLE_OVERFLOW) {
			if self.0<=0 {
				panic!("PseudoDouble Log2 of non-positive number");
			}
		}
		let exponent=self.0&EXP_MASK;
		let e=exponent-PSEUDO_DOUBLE_EXP_BIAS-2;
		let mantissa=((self.0&EXP_MASK_INV)<<2)>>1;
		let log_frac=log2_64_fixed(mantissa as u64);
		if e==0 {
			if log_frac==0 {
				return PseudoDouble(0);
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
	//return (((e<<(PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65))+(log_frac>>(64-lead_bits)))&EXP_MASK_INV)+PSEUDO_DOUBLE_EXP_BIAS+65-lead_bits;
	}
}
/*
	inline PseudoDouble(double f) {val=double_to_pdi(f);}
	inline operator double() const {return pdi_to_double(val);}
	static PseudoDouble create(pseudo_double_i pdi) {PseudoDouble ret;ret.val=pdi;return ret;}
	friend PseudoDouble exp(const PseudoDouble x);
	friend PseudoDouble log2(const PseudoDouble x);
	friend PseudoDouble log(const PseudoDouble x);
	friend PseudoDouble log10(const PseudoDouble x);
	friend PseudoDouble pow(const PseudoDouble x, const PseudoDouble y);
	friend PseudoDouble sin_rev(const PseudoDouble x);
	friend PseudoDouble cos_rev(const PseudoDouble x);
	friend PseudoDouble atan2_rev(const PseudoDouble y, const PseudoDouble x);
	friend PseudoDouble sin(const PseudoDouble x);
	friend PseudoDouble cos(const PseudoDouble x);
	friend PseudoDouble atan2(const PseudoDouble y, const PseudoDouble x);
	friend PseudoDouble abs(const PseudoDouble x);
	friend PseudoDouble PD_create_fixed10(int64_t x, int32_t e);
	friend PseudoDouble PD_create_fixed2(int64_t x, int32_t e);
	friend int64_t PD_get_fixed2(PseudoDouble x, int32_t e);
	friend PseudoDouble PD_from_string(std::string str);
*/

fn main() {
	let x=PseudoDouble::from(2);
	let y=PseudoDouble::from(7);
	let z=i64::from(x+y);
	println!("z={}",z);
}
