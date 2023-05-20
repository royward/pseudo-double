# Functions

The form of these is:

### \<function name\>

	mlib: <standard math library signature>
	C: <C signature>
	C++: <C++ signature>
	<optional description>

## Basic Arithmetic

### Equal

	mlib: x==y
	C: (use integer x==y)
	inline bool PseudoFloat::operator==(const PseudoFloat x) const;

### Not equal

	mlib: x!=y
	C: (use integer x!=y)
	inline bool PseudoFloat::operator!=(const PseudoFloat x) const;

### Greater than or equal to

	mlib: x>=y
	C: inline int pf_gte(pseudo_float x, pseudo_float y);
	inline bool PseudoFloat::operator>=(const PseudoFloat x) const;

### Less than

	mlib: x<y
	C: (use pf_gt with reversed arguments)
	inline bool PseudoFloat::operator<(const PseudoFloat x) const;

### Less than or equal to

	mlib: x<=y
	C: (use pf_gte with reversed arguments)
	inline bool PseudoFloat::operator<=(const PseudoFloat x) const;

## Quick comparison with 0. For C, just do a straight integer comparison with zero

	C++: inline bool PseudoFloat::gt_zero() const;
	C++: inline bool PseudoFloat::gte_zero() const;
	C++: inline bool PseudoFloat::lt_zero() const;
	C++: inline bool PseudoFloat::lte_zero() const;
	C++: inline bool PseudoFloat::eq_zero() const;
	C++: inline bool PseudoFloat::neq_zero() const;

## Conversion

### convert double to pseudo-float

	C: pseudo_float double_to_pf(double d);
	C++: inline PseudoFloat(double f);

### convert signed integer to pseudo-float

	C: pseudo_float int64_to_pf(int64_t d);
	C++: inline PseudoFloat(int16_t f);
	C++: inline PseudoFloat(int32_t f);
	C++: inline PseudoFloat(int64_t f);

### convert unsigned integer to pseudo-float

	C: pseudo_float uint64_to_pf(uint64_t d);
	C++: inline PseudoFloat(uint16_t f);
	C++: inline PseudoFloat(uint32_t f);
	C++: inline PseudoFloat(uint64_t f);

### convert pseudo-float to double

	C: double pf_to_double(pseudo_float d);
	C++: inline operator double() const {return pf_to_double(val);}

### convert pseudo-float to signed integer

	C: int64_t pf_to_int64(pseudo_float d);
	C++: inline operator PseudoFloat::int16_t() const;
	C++: inline operator PseudoFloat::int32_t() const;
	C++: inline operator PseudoFloat::int64_t() const;

### convert pseudo-float to unsigned integer

	C: uint64_t pf_to_uint64(pseudo_float d);
	C++: inline operator PseudoFloat::uint16_t() const;
	C++: inline operator PseudoFloat::uint32_t() const;
	C++: inline operator PseudoFloat::uint64_t() const;

	C: pseudo_float sint64fixed2_to_pf(int64_t d, int32_t e);
	C++: int64_t pf_to_sint64fixed2(pseudo_float d, int32_t e);
	C++: inline pseudo_float sint64fixed10_to_pf(int64_t d, int32_t e);

## Miscellaneous

### floor

	mlib: floor(x)
	C: pseudo_float pf_floor(pseudo_float x);
	C++: PseudoFloat floor(const PseudoFloat x);

### ceiling

	mlib: ceil(x)
	C: pseudo_float pf_ceil(pseudo_float x);
	C++: PseudoFloat ceil(const PseudoFloat x);

### round

	mlib: round(x)
	C: pseudo_float pf_round(pseudo_float x);
	C++: PseudoFloat round(const PseudoFloat x);
	Returns the integral value that is nearest to x, with halfway cases rounded away from zero.

### ldexp

	mlib: ldexp(x)
	C: inline pseudo_float pf_ldexp(pseudo_float x, int y);
	C++: PseudoFloat ldexp(const PseudoFloat x, int y);

## Square root

### square root

	mlib: sqrt(x)
	C: pseudo_float pf_sqrt(pseudo_float x);
	C++: PseudoFloat sqrt(const PseudoFloat x);

### inverse square root

	mlib: 1.0/sqrt(x)
	C: pseudo_float pf_inv_sqrt(pseudo_float x);
	C++: PseudoFloat inv_sqrt(const PseudoFloat x);
	This is faster than sqrt and useful for normalizing

## Helper functions

	PseudoFloat PF_create_fixed10(int64_t x, int32_t e);
	PseudoFloat PF_create_fixed2(int64_t x, int32_t e);
	int64_t PF_get_fixed2(PseudoFloat x, int32_t e);

	C++: inline pseudo_float PseudoFloat::get_internal() const;

	C++: inline void PseudoFloat::set_internal(pseudo_float f);

## Exponential and power functions

The base fuctions are pf_exp2 and pf_log2, so:

* pf_exp2(n) is exactly 2^n (subject to no overflow)
* pf_log2(2^n) is exactly n

### binary exponent

	mlib: exp2(x)
	C: pseudo_float pf_exp2(pseudo_float x);
	C++: PseudoFloat exp2(const PseudoFloat x);

### logarithm base 2

	mlib: log2(x)
	C: pseudo_float pf_log2(pseudo_float x);
	C++: PseudoFloat log2(const PseudoFloat x);

### exponential function

	mlib: exp(x)
	C: pseudo_float pf_exp(pseudo_float x);
	C++: PseudoFloat exp(const PseudoFloat x);

### logarithm base e

	mlib: log(x)
	C: pseudo_float pf_log(pseudo_float x);
	C++: PseudoFloat log(const PseudoFloat x);

### logarithm base 10

	mlib: log10(x)
	C: pseudo_float pf_log10(pseudo_float x);
	C++: PseudoFloat log10(const PseudoFloat x);

### power, x^y

	mlib: pow(x,y)
	C: pseudo_float pf_pow(pseudo_float x, pseudo_float y);
	C++: PseudoFloat pow(const PseudoFloat x, const PseudoFloat y);

## Trigonometry

Becasuse this is aimed at games rather than scientific applications, the base trigonometry functions use revolutions (turns). A revolution is 2$\pi$ radians.

The following identities are exactly true for all integer n and pseudofloat x>0:

* pf_sin_rev(n/2)=0
* pf_sin_rev(n+1/4)=1
* pf_sin_rev(n+3/4)=-1
* pf_cos_rev(n/2+1/4)=0
* pf_cos_rev(n)=1
* pf_cos_rev(n+1/2)=-1
* pf_atan2_rev(0,0)  =0
* pf_atan2_rev(0,x)  =0
* pf_atan2_rev(x,x)  =1/8
* pf_atan2_rev(x,0)  =1/4
* pf_atan2_rev(x,-x) =3/8
* pf_atan2_rev(0,-x) =1/2
* pf_atan2_rev(-x,-x)=5/8
* pf_atan2_rev(-x,0) =3/4
* pf_atan2_rev(-x,x) =7/8

### sin of revolutions

	mlib: sin(x*(2.0*M_PI))
	C: pseudo_float pf_sin_rev(pseudo_float x);
	C++: PseudoFloat sin_rev(const PseudoFloat x);

### cos of revolutions

	mlib: cos(x*(2.0*M_PI))
	C: pseudo_float pf_cos_rev(pseudo_float x);
	C++: PseudoFloat cos_rev(const PseudoFloat x);

### atan2 in revolutions

	mlib: atan2(y,x)/(2.0*M_PI)
	C: pseudo_float pf_atan2_rev(pseudo_float y, pseudo_float x);
	C++: PseudoFloat atan2_rev(const PseudoFloat y, const PseudoFloat x);

### sin of radians

	mlib: sin(x)
	C: pseudo_float pf_sin(pseudo_float x);
	C++: PseudoFloat sin(const PseudoFloat x);

### cos of radians

	mlib: cos(x)
	C: pseudo_float pf_cos(pseudo_float x);
	C++: PseudoFloat cos(const PseudoFloat x);

### atan2 in radians

	mlib: atan2(y,x)
	C: pseudo_float pf_atan2(pseudo_float y, pseudo_float x);
	C++: PseudoFloat atan2(const PseudoFloat y, const PseudoFloat x);

## Fixed integer helper functions

	C/C++: inline uint64_t multu64hi(uint64_t x,uint64_t y);
	calculate (x*y)>>64
	Useful if x,y and result are considered 0.64 fixed ints

	C/C++: uint64_t inv_sqrt64_fixed(uint64_t x);
	calculate 1/sqrt(x)
	x is a 2.62 unsigned fixed in the range (1,4)
	result is 1.63 unsigned fixed in the range (0.5,1)

	C/C++: uint64_t exp2_64_fixed(uint64_t x);
	calculate 2^x
	x is a 0.64 unsigned fixed in the range [0,1)
	result is 2.62 unsigned fixed in the range [1,2)

	C/C++: uint64_t log2_64_fixed(uint64_t x);
	calculate ln2(x+1)
	x is a 1.63 unsigned fixed in the range [0,1)
	result is 1.63 unsigned fixed in the range [0,1)

	C/C++: uint64_t sin_rev_64_fixed(uint64_t x);
	calculate sin(x)
	x is a 0.64 unsigned fixed in the range [0,1/4] in revolutions
	result is 2.62 unsigned fixed in the range [0,1]

	C/C++: uint64_t atan_rev_64_fixed(uint64_t x);
	calculate atan(x)
	x is a 2.62 unsigned fixed in the range [0,1]
	result is (-1).65 unsigned fixed in the range [0,1/8] in revolutions
	(-1).65 means the values 0x0000000000000000..0xFFFFFFFFFFFFFFFF maps to 0..1/2
