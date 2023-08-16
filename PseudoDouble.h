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

#pragma once
#ifndef PSEUDO_DOUBLE_CPP_H
#define PSEUDO_DOUBLE_CPP_H

#include <stdexcept>

#ifndef PD_DO_ERROR_OVERFLOW
#define PD_DO_ERROR_OVERFLOW throw std::overflow_error("overflow");
#endif

#ifndef PD_DO_ERROR_RANGE
#define PD_DO_ERROR_RANGE throw std::range_error("range");
#endif

#include "pseudo_double.h"

class PseudoDouble {
public:
	PseudoDouble() {val=0;};
	inline PseudoDouble(double f) {val=double_to_pdi(f);}
	inline PseudoDouble(int16_t f) {val=int64_to_pdi(f);}
	inline PseudoDouble(int32_t f) {val=int64_to_pdi(f);}
	inline PseudoDouble(int64_t f) {val=int64_to_pdi(f);}
	inline PseudoDouble(uint16_t f) {val=uint64_to_pdi(f);}
	inline PseudoDouble(uint32_t f) {val=uint64_to_pdi(f);}
	inline PseudoDouble(uint64_t f) {val=uint64_to_pdi(f);}
#ifndef _MSC_VER // windows
	inline PseudoDouble(long long f) { val = int64_to_pdi(f); }
	inline PseudoDouble(unsigned long long f) {val=uint64_to_pdi(f);}
#endif
	inline operator double() const {return pdi_to_double(val);}
	inline operator int16_t() const {return static_cast<int16_t>(pdi_to_int64(val));}
	inline operator int32_t() const {return static_cast<int32_t>(pdi_to_int64(val));}
	inline operator int64_t() const {return pdi_to_int64(val);}
	inline operator uint16_t() const {return static_cast<uint16_t>(pdi_to_uint64(val));}
	inline operator uint32_t() const {return static_cast<uint32_t>(pdi_to_uint64(val));}
	inline operator uint64_t() const {return pdi_to_uint64(val);}
	inline PseudoDouble operator-() const {return PseudoDouble::create(pdi_neg(val));}
	inline PseudoDouble operator+(const PseudoDouble x) const {return PseudoDouble::create(pdi_add(val,x.val));}
	inline PseudoDouble operator-(const PseudoDouble x) const {return PseudoDouble::create(pdi_sub(val,x.val));}
	inline PseudoDouble operator*(const PseudoDouble x) const {return PseudoDouble::create(pdi_mult(val,x.val));}
	inline PseudoDouble operator/(const PseudoDouble x) const {return PseudoDouble::create(pdi_div(val,x.val));}
	inline PseudoDouble operator+=(const PseudoDouble x) {val=pdi_add(val,x.val); return *this;}
	inline PseudoDouble operator-=(const PseudoDouble x) {val=pdi_sub(val,x.val); return *this;}
	inline PseudoDouble operator*=(const PseudoDouble x) {val=pdi_mult(val,x.val); return *this;}
	inline PseudoDouble operator/=(const PseudoDouble x) {val=pdi_div(val,x.val); return *this;}
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
	inline pseudo_double_i get_internal() const {return val;}
	inline void set_internal(pseudo_double_i f) {val=f;}
private:
	static PseudoDouble create(pseudo_double_i pdi) {PseudoDouble ret;ret.val=pdi;return ret;}
	pseudo_double_i val;
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
};

inline PseudoDouble floor(const PseudoDouble x) {return PseudoDouble::create(pdi_floor(x.val));}
inline PseudoDouble ceil(const PseudoDouble x) {return PseudoDouble::create(pdi_ceil(x.val));}
inline PseudoDouble round(const PseudoDouble x) {return PseudoDouble::create(pdi_round(x.val));}
inline PseudoDouble sqrt(const PseudoDouble x) {return PseudoDouble::create(pdi_sqrt(x.val));}
inline PseudoDouble inv_sqrt(const PseudoDouble x) {return PseudoDouble::create(pdi_inv_sqrt(x.val));}
inline PseudoDouble ldexp(const PseudoDouble x, int y) {return PseudoDouble::create(pdi_ldexp(x.val,y));}
inline PseudoDouble exp2(const PseudoDouble x) {return PseudoDouble::create(pdi_exp2(x.val));}
inline PseudoDouble exp(const PseudoDouble x) {return PseudoDouble::create(pdi_exp(x.val));}
inline PseudoDouble log2(const PseudoDouble x) {return PseudoDouble::create(pdi_log2(x.val));}
inline PseudoDouble log(const PseudoDouble x) {return PseudoDouble::create(pdi_log(x.val));}
inline PseudoDouble log10(const PseudoDouble x) {return PseudoDouble::create(pdi_log10(x.val));}
inline PseudoDouble pow(const PseudoDouble x, const PseudoDouble y) {return PseudoDouble::create(pdi_pow(x.val,y.val));}
inline PseudoDouble sin_rev(const PseudoDouble x) {return PseudoDouble::create(pdi_sin_rev(x.val));}
inline PseudoDouble cos_rev(const PseudoDouble x) {return PseudoDouble::create(pdi_cos_rev(x.val));}
inline PseudoDouble atan2_rev(const PseudoDouble y, const PseudoDouble x) {return PseudoDouble::create(pdi_atan2_rev(y.val,x.val));}
inline PseudoDouble sin(const PseudoDouble x) {return PseudoDouble::create(pdi_sin(x.val));}
inline PseudoDouble cos(const PseudoDouble x) {return PseudoDouble::create(pdi_cos(x.val));}
inline PseudoDouble atan2(const PseudoDouble y, const PseudoDouble x) {return PseudoDouble::create(pdi_atan2(y.val,x.val));}
inline PseudoDouble abs(const PseudoDouble x) {return PseudoDouble::create(pdi_abs(x.val));}
inline PseudoDouble fabs(const PseudoDouble x) {return PseudoDouble::create(pdi_abs(x.val));}
inline PseudoDouble PD_create_fixed10(int64_t x, int32_t e) {return PseudoDouble::create(int64fixed10_to_pdi(x,e));}
inline PseudoDouble PD_create_fixed2(int64_t x, int32_t e) {return PseudoDouble::create(int64fixed2_to_pdi(x,e));}
inline int64_t PD_get_fixed2(PseudoDouble x, int32_t e) {return pdi_to_int64fixed2(x.val,e);}

const static PseudoDouble PD_HALF=PD_create_fixed2(1,-1);
const static PseudoDouble PD_ZERO=PseudoDouble(0);
const static PseudoDouble PD_ONE=PseudoDouble(1U);
const static PseudoDouble PD_TWO=PseudoDouble(2U);
const static PseudoDouble PD_PI_DIV_2=PD_create_fixed10(1570796326794896619,-18);
const static PseudoDouble PD_PI=PD_create_fixed10(3141592653589793238,-18);
const static PseudoDouble PD_TAU=PD_create_fixed10(6283185307179586477,-18);
const static PseudoDouble PD_2_DIV_PI=PD_ONE/PD_PI_DIV_2;
const static PseudoDouble PD_INV_PI=PD_ONE/PD_PI;
const static PseudoDouble PD_INV_TAU=PD_ONE/PD_TAU;

#endif // PSEUDO_DOUBLE_CPP_H
