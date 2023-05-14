#pragma once
#ifndef PSEUDO_FLOAT_CPP_H
#define PSEUDO_FLOAT_CPP_H

#include <stdexcept>

#ifndef PF_DO_ERROR_OVERFLOW
#define PF_DO_ERROR_OVERFLOW throw std::overflow_error("overflow");
#endif

#ifndef PF_DO_ERROR_RANGE
#define PF_DO_ERROR_RANGE throw std::range_error("range");
#endif

#include "pseudo_float.h"

class PseudoFloat {
public:
	PseudoFloat() {val=0;};
	inline PseudoFloat(double f) {val=double_to_pf(f);}
	inline PseudoFloat(int16_t f) {val=int64_to_pf(f);}
	inline PseudoFloat(int32_t f) {val=int64_to_pf(f);}
	inline PseudoFloat(int64_t f) {val=int64_to_pf(f);}
	inline PseudoFloat(uint16_t f) {val=uint64_to_pf(f);}
	inline PseudoFloat(uint32_t f) {val=uint64_to_pf(f);}
	inline PseudoFloat(uint64_t f) {val=uint64_to_pf(f);}
	inline operator double() const {return pf_to_double(val);}
	inline operator int16_t() const {return pf_to_int64(val);}
	inline operator int32_t() const {return pf_to_int64(val);}
	inline operator int64_t() const {return pf_to_int64(val);}
	inline operator uint16_t() const {return pf_to_uint64(val);}
	inline operator uint32_t() const {return pf_to_uint64(val);}
	inline operator uint64_t() const {return pf_to_uint64(val);}
	inline PseudoFloat operator-() const {return PseudoFloat::create(pf_neg(val));}
	inline PseudoFloat operator+(const PseudoFloat x) const {return PseudoFloat::create(pf_add(val,x.val));}
	inline PseudoFloat operator-(const PseudoFloat x) const {return PseudoFloat::create(pf_sub(val,x.val));}
	inline PseudoFloat operator*(const PseudoFloat x) const {return PseudoFloat::create(pf_mult(val,x.val));}
	inline PseudoFloat operator/(const PseudoFloat x) const {return PseudoFloat::create(pf_div(val,x.val));}
	inline PseudoFloat operator+=(const PseudoFloat x) {val=pf_add(val,x.val); return *this;}
	inline PseudoFloat operator-=(const PseudoFloat x) {val=pf_sub(val,x.val); return *this;}
	inline PseudoFloat operator*=(const PseudoFloat x) {val=pf_mult(val,x.val); return *this;}
	inline PseudoFloat operator/=(const PseudoFloat x) {val=pf_div(val,x.val); return *this;}
	inline bool operator==(const PseudoFloat x) const {return val==x.val;}
	inline bool operator!=(const PseudoFloat x) const {return val!=x.val;}
	inline bool operator>(const PseudoFloat x) const {return pf_gt(val,x.val);}
	inline bool operator>=(const PseudoFloat x) const {return pf_gte(val,x.val);}
	inline bool operator<(const PseudoFloat x) const {return pf_gt(x.val,val);}
	inline bool operator<=(const PseudoFloat x) const {return pf_gte(x.val,val);}
	inline bool gt_zero() const {return ((signed_pf_internal)val)>(signed_pf_internal)0;}
	inline bool gte_zero() const {return ((signed_pf_internal)val)>=(signed_pf_internal)0;}
	inline bool lt_zero() const {return ((signed_pf_internal)val)<(signed_pf_internal)0;}
	inline bool lte_zero() const {return ((signed_pf_internal)val)<=(signed_pf_internal)0;}
	inline bool eq_zero() const {return ((signed_pf_internal)val)==(signed_pf_internal)0;}
	inline bool neq_zero() const {return ((signed_pf_internal)val)!=(signed_pf_internal)0;}
	inline pseudo_float get_internal() const {return val;}
	inline void set_internal(pseudo_float f) {val=f;}
private:
	static PseudoFloat create(pseudo_float pf) {PseudoFloat ret;ret.val=pf;return ret;}
	pseudo_float val;
	friend PseudoFloat floor(const PseudoFloat x);
	friend PseudoFloat ceil(const PseudoFloat x);
	friend PseudoFloat round(const PseudoFloat x);
	friend PseudoFloat sqrt(const PseudoFloat x);
	friend PseudoFloat inv_sqrt(const PseudoFloat x);
	friend PseudoFloat ldexp(const PseudoFloat x, int y);
	friend PseudoFloat exp2(const PseudoFloat x);
	friend PseudoFloat exp(const PseudoFloat x);
	friend PseudoFloat log2(const PseudoFloat x);
	friend PseudoFloat log(const PseudoFloat x);
	friend PseudoFloat log10(const PseudoFloat x);
	friend PseudoFloat pow(const PseudoFloat x, const PseudoFloat y);
	friend PseudoFloat sin_rev(const PseudoFloat x);
	friend PseudoFloat cos_rev(const PseudoFloat x);
	friend PseudoFloat atan2_rev(const PseudoFloat y, const PseudoFloat x);
	friend PseudoFloat sin(const PseudoFloat x);
	friend PseudoFloat cos(const PseudoFloat x);
	friend PseudoFloat atan2(const PseudoFloat y, const PseudoFloat x);
	friend PseudoFloat abs(const PseudoFloat x);
	friend PseudoFloat fabs(const PseudoFloat x);
	friend PseudoFloat PF_create_fixed10(int64_t x, int32_t e);
	friend PseudoFloat PF_create_fixed2(int64_t x, int32_t e);
	friend int64_t PF_get_fixed2(PseudoFloat x, int32_t e);
};

inline PseudoFloat floor(const PseudoFloat x) {return PseudoFloat::create(pf_floor(x.val));}
inline PseudoFloat ceil(const PseudoFloat x) {return PseudoFloat::create(pf_ceil(x.val));}
inline PseudoFloat round(const PseudoFloat x) {return PseudoFloat::create(pf_round(x.val));}
inline PseudoFloat sqrt(const PseudoFloat x) {return PseudoFloat::create(pf_sqrt(x.val));}
inline PseudoFloat inv_sqrt(const PseudoFloat x) {return PseudoFloat::create(pf_inv_sqrt(x.val));}
inline PseudoFloat ldexp(const PseudoFloat x, int y) {return PseudoFloat::create(pf_ldexp(x.val,y));}
inline PseudoFloat exp2(const PseudoFloat x) {return PseudoFloat::create(pf_exp2(x.val));}
inline PseudoFloat exp(const PseudoFloat x) {return PseudoFloat::create(pf_exp(x.val));}
inline PseudoFloat log2(const PseudoFloat x) {return PseudoFloat::create(pf_log2(x.val));}
inline PseudoFloat log(const PseudoFloat x) {return PseudoFloat::create(pf_log(x.val));}
inline PseudoFloat log10(const PseudoFloat x) {return PseudoFloat::create(pf_log10(x.val));}
inline PseudoFloat pow(const PseudoFloat x, const PseudoFloat y) {return PseudoFloat::create(pf_pow(x.val,y.val));}
inline PseudoFloat sin_rev(const PseudoFloat x) {return PseudoFloat::create(pf_sin_rev(x.val));}
inline PseudoFloat cos_rev(const PseudoFloat x) {return PseudoFloat::create(pf_cos_rev(x.val));}
inline PseudoFloat atan2_rev(const PseudoFloat y, const PseudoFloat x) {return PseudoFloat::create(pf_atan2_rev(y.val,x.val));}
inline PseudoFloat sin(const PseudoFloat x) {return PseudoFloat::create(pf_sin(x.val));}
inline PseudoFloat cos(const PseudoFloat x) {return PseudoFloat::create(pf_cos(x.val));}
inline PseudoFloat atan2(const PseudoFloat y, const PseudoFloat x) {return PseudoFloat::create(pf_atan2(y.val,x.val));}
inline PseudoFloat abs(const PseudoFloat x) {return PseudoFloat::create(pf_abs(x.val));}
inline PseudoFloat fabs(const PseudoFloat x) {return PseudoFloat::create(pf_abs(x.val));}
inline PseudoFloat PF_create_fixed10(int64_t x, int32_t e) {return PseudoFloat::create(int64fixed10_to_pf(x,e));}
inline PseudoFloat PF_create_fixed2(int64_t x, int32_t e) {return PseudoFloat::create(int64fixed2_to_pf(x,e));}
inline int64_t PF_get_fixed2(PseudoFloat x, int32_t e) {return pf_to_int64fixed2(x.val,e);}

const static PseudoFloat PF_HALF=PF_create_fixed2(1,-1);
const static PseudoFloat PF_ZERO=PseudoFloat(0);
const static PseudoFloat PF_ONE=PseudoFloat(1U);
const static PseudoFloat PF_TWO=PseudoFloat(2U);
const static PseudoFloat PF_PI_DIV_2=PF_create_fixed10(1570796326794896619,-18);
const static PseudoFloat PF_PI=PF_create_fixed10(3141592653589793238,-18);
const static PseudoFloat PF_TAU=PF_create_fixed10(6283185307179586477,-18);
const static PseudoFloat PF_2_DIV_PI=PF_ONE/PF_PI_DIV_2;
const static PseudoFloat PF_INV_PI=PF_ONE/PF_PI;
const static PseudoFloat PF_INV_TAU=PF_ONE/PF_TAU;

#endif
