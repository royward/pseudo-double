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
use std::convert::From;

#[derive(Debug, Copy, Clone)]
struct PseudoDouble(i64);

const PSEUDO_DOUBLE_TOTAL_BITS: i32 = 64;
const PSEUDO_DOUBLE_EXP_BITS: i32 = 16;
const EXP_MASK: i64 = (1<<PSEUDO_DOUBLE_EXP_BITS)-1;
const EXP_MASK_INV: i64 = !EXP_MASK as i64;
const PSEUDO_DOUBLE_HALF_ULP: i64 = (1<<(PSEUDO_DOUBLE_EXP_BITS-1))-1;
const PSEUDO_DOUBLE_EXP_BIAS: i64 = 1<<(PSEUDO_DOUBLE_EXP_BITS-1);

fn shift_left_signed(x:i64, shift:i32) -> i64 {
	if shift>0 {x<<shift} else {x>>-shift}
}

fn shift_left_unsigned(x:u64, shift:i32) -> u64 {
	if shift>0 {x<<shift} else {x>>-shift}
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

/*
	inline PseudoDouble(double f) {val=double_to_pdi(f);}
	inline operator double() const {return pdi_to_double(val);}
	inline bool operator==(const PseudoDouble x) const {return val==x.val;}
	inline bool operator!=(const PseudoDouble x) const {return val!=x.val;}
	inline bool operator>(const PseudoDouble x) const {return pdi_gt(val,x.val);}
	inline bool operator>=(const PseudoDouble x) const {return pdi_gte(val,x.val);}
	inline bool operator<(const PseudoDouble x) const {return pdi_gt(x.val,val);}
	inline bool operator<=(const PseudoDouble x) const {return pdi_gte(x.val,val);}
	inline bool gt_zero() const {return ((signed_pd_internal)val)>(signed_pd_internal)0;}
	inline bool gte_zero() const {return ((signed_pd_internal)val)>=(signed_pd_internal)0;}
	inline bool lt_zero() const {return ((signed_pd_internal)val)<(signed_pd_internal)0;}
	inline bool lte_zero() const {return ((signed_pd_internal)val)<=(signed_pd_internal)0;}
	inline bool eq_zero() const {return ((signed_pd_internal)val)==(signed_pd_internal)0;}
	inline bool neq_zero() const {return ((signed_pd_internal)val)!=(signed_pd_internal)0;}
	static PseudoDouble create(pseudo_double_i pdi) {PseudoDouble ret;ret.val=pdi;return ret;}
	friend PseudoDouble min(const PseudoDouble x, const PseudoDouble y);
	friend PseudoDouble max(const PseudoDouble x, const PseudoDouble y);
	friend PseudoDouble floor(const PseudoDouble x);
	friend PseudoDouble ceil(const PseudoDouble x);
	friend PseudoDouble round(const PseudoDouble x);
	friend PseudoDouble sqrt(const PseudoDouble x);
	friend PseudoDouble inv_sqrt(const PseudoDouble x);
	friend PseudoDouble ldexp(const PseudoDouble x, int y);
	friend PseudoDouble exp2(const PseudoDouble x);
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
	friend PseudoDouble fabs(const PseudoDouble x);
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
