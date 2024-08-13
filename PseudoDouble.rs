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

use std::ops::{Add, AddAssign};
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

impl From<i32> for PseudoDouble {
    fn from(x : i32) -> Self {
		return PseudoDouble::from(x as i64);
    }
}

impl From<u64> for PseudoDouble {
    fn from(x : u64) -> Self {
		if x==0 {
			return PseudoDouble(0);
		} else {
			let lead_bits=x.leading_zeros() as i32;
			return PseudoDouble(((shift_left_unsigned(x,PSEUDO_DOUBLE_TOTAL_BITS+lead_bits-65)) as i64&EXP_MASK_INV) +PSEUDO_DOUBLE_EXP_BIAS as i64+65-lead_bits as i64);
		}
    }
}

impl From<u32> for PseudoDouble {
    fn from(x : u32) -> Self {
		return PseudoDouble::from(x as u64);
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

impl From<PseudoDouble> for u32 {
    fn from(x: PseudoDouble) -> Self {
		return u64::from(x) as u32;
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

impl From<PseudoDouble> for i32 {
    fn from(x: PseudoDouble) -> Self {
		return i64::from(x) as i32;
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
			if new_exponent<0 {
				return PseudoDouble(0);
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

impl AddAssign for PseudoDouble {
    fn add_assign(&mut self, other: Self) {
        *self = *self + other;
    }
}

fn main() {
	let x=PseudoDouble::from(2);
	let y=PseudoDouble::from(7);
	let z=i64::from(x+y);
	println!("z={}",z);
}
