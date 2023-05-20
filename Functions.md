# Functions

The form of these is:

### \<function name\>

	mlib: <standard math library or built in signature>
	C: <C signature>
	C++: <C++ signature>
	<optional description>

## Basic Arithmetic

### Add

	mlib: x+y
	C: inline int pd_add(pseudo_double x, pseudo_double y);
	inline PseudoDouble PseudoDouble::operator+(const PseudoDouble x) const;

### Subtract

	mlib: x-y
	C: inline int pd_sub(pseudo_double x, pseudo_double y);
	inline PseudoDouble PseudoDouble::operator-(const PseudoDouble x) const;

### Multiply

	mlib: x*y
	C: inline int pd_mult(pseudo_double x, pseudo_double y);
	inline PseudoDouble PseudoDouble::operator*(const PseudoDouble x) const;

### Divide

	mlib: x/y
	C: inline int pd_div(pseudo_double x, pseudo_double y);
	inline PseudoDouble PseudoDouble::operator/(const PseudoDouble x) const;

### Negate

	mlib: -x
	C: inline int pd_neg(pseudo_double x);
	inline PseudoDouble PseudoDouble::operator-() const;

### Equal

	mlib: x==y
	C: (use integer x==y)
	inline bool PseudoDouble::operator==(const PseudoDouble x) const;

### Not equal

	mlib: x!=y
	C: (use integer x!=y)
	inline bool PseudoDouble::operator!=(const PseudoDouble x) const;

### Greater than or equal to

	mlib: x>=y
	C: inline int pd_gte(pseudo_double x, pseudo_double y);
	inline bool PseudoDouble::operator>=(const PseudoDouble x) const;

### Less than

	mlib: x<y
	C: (use pd_gt with reversed arguments)
	inline bool PseudoDouble::operator<(const PseudoDouble x) const;

### Less than or equal to

	mlib: x<=y
	C: (use pd_gte with reversed arguments)
	inline bool PseudoDouble::operator<=(const PseudoDouble x) const;

## Quick comparison with 0.

For C, just do a straight integer comparison with zero

	C++: inline bool PseudoDouble::gt_zero() const;
	C++: inline bool PseudoDouble::gte_zero() const;
	C++: inline bool PseudoDouble::lt_zero() const;
	C++: inline bool PseudoDouble::lte_zero() const;
	C++: inline bool PseudoDouble::eq_zero() const;
	C++: inline bool PseudoDouble::neq_zero() const;

## Conversion

### convert double to pseudo-double

	C: pseudo_double double_to_pd(double d);
	C++: inline PseudoDouble(double f);

### convert signed integer to pseudo-double

	C: pseudo_double int64_to_pd(int64_t d);
	C++: inline PseudoDouble(int16_t f);
	C++: inline PseudoDouble(int32_t f);
	C++: inline PseudoDouble(int64_t f);

### convert unsigned integer to pseudo-double

	C: pseudo_double uint64_to_pd(uint64_t d);
	C++: inline PseudoDouble(uint16_t f);
	C++: inline PseudoDouble(uint32_t f);
	C++: inline PseudoDouble(uint64_t f);

### convert pseudo-double to double

	C: double pd_to_double(pseudo_double d);
	C++: inline operator double() const {return pd_to_double(val);}

### convert pseudo-double to signed integer

	C: int64_t pd_to_int64(pseudo_double d);
	C++: inline operator PseudoDouble::int16_t() const;
	C++: inline operator PseudoDouble::int32_t() const;
	C++: inline operator PseudoDouble::int64_t() const;

### convert pseudo-double to unsigned integer

	C: uint64_t pd_to_uint64(pseudo_double d);
	C++: inline operator PseudoDouble::uint16_t() const;
	C++: inline operator PseudoDouble::uint32_t() const;
	C++: inline operator PseudoDouble::uint64_t() const;

	C: pseudo_double sint64fixed2_to_pd(int64_t d, int32_t e);
	C++: int64_t pd_to_sint64fixed2(pseudo_double d, int32_t e);
	C++: inline pseudo_double sint64fixed10_to_pd(int64_t d, int32_t e);

## Miscellaneous

### floor

	mlib: floor(x)
	C: pseudo_double pd_floor(pseudo_double x);
	C++: PseudoDouble floor(const PseudoDouble x);

### ceiling

	mlib: ceil(x)
	C: pseudo_double pd_ceil(pseudo_double x);
	C++: PseudoDouble ceil(const PseudoDouble x);

### round

	mlib: round(x)
	C: pseudo_double pd_round(pseudo_double x);
	C++: PseudoDouble round(const PseudoDouble x);
	Returns the integral value that is nearest to x, with halfway cases rounded away from zero.

### ldexp

	mlib: ldexp(x)
	C: inline pseudo_double pd_ldexp(pseudo_double x, int y);
	C++: PseudoDouble ldexp(const PseudoDouble x, int y);

## Square root

### square root

	mlib: sqrt(x)
	C: pseudo_double pd_sqrt(pseudo_double x);
	C++: PseudoDouble sqrt(const PseudoDouble x);

### inverse square root

	mlib: 1.0/sqrt(x)
	C: pseudo_double pd_inv_sqrt(pseudo_double x);
	C++: PseudoDouble inv_sqrt(const PseudoDouble x);
	This is faster than sqrt and useful for normalizing

## Helper functions

	PseudoDouble PD_create_fixed10(int64_t x, int32_t e);
	PseudoDouble PD_create_fixed2(int64_t x, int32_t e);
	int64_t PD_get_fixed2(PseudoDouble x, int32_t e);

	C++: inline pseudo_double PseudoDouble::get_internal() const;

	C++: inline void PseudoDouble::set_internal(pseudo_double f);

## Exponential and power functions

The base fuctions are pd_exp2 and pd_log2, so:

* pd_exp2(n) is exactly 2^n (subject to no overflow)
* pd_log2(2^n) is exactly n

### binary exponent

	mlib: exp2(x)
	C: pseudo_double pd_exp2(pseudo_double x);
	C++: PseudoDouble exp2(const PseudoDouble x);

### logarithm base 2

	mlib: log2(x)
	C: pseudo_double pd_log2(pseudo_double x);
	C++: PseudoDouble log2(const PseudoDouble x);

### exponential function

	mlib: exp(x)
	C: pseudo_double pd_exp(pseudo_double x);
	C++: PseudoDouble exp(const PseudoDouble x);

### logarithm base e

	mlib: log(x)
	C: pseudo_double pd_log(pseudo_double x);
	C++: PseudoDouble log(const PseudoDouble x);

### logarithm base 10

	mlib: log10(x)
	C: pseudo_double pd_log10(pseudo_double x);
	C++: PseudoDouble log10(const PseudoDouble x);

### power, x^y

	mlib: pow(x,y)
	C: pseudo_double pd_pow(pseudo_double x, pseudo_double y);
	C++: PseudoDouble pow(const PseudoDouble x, const PseudoDouble y);

## Trigonometry

Becasuse this is aimed at games rather than scientific applications, the base trigonometry functions use revolutions (turns). A revolution is 2$\pi$ radians.

The following identities are exactly true for all integer n and pseudodouble x>0:

* pd_sin_rev(n/2)=0
* pd_sin_rev(n+1/4)=1
* pd_sin_rev(n+3/4)=-1
* pd_cos_rev(n/2+1/4)=0
* pd_cos_rev(n)=1
* pd_cos_rev(n+1/2)=-1
* pd_atan2_rev(0,0)  =0
* pd_atan2_rev(0,x)  =0
* pd_atan2_rev(x,x)  =1/8
* pd_atan2_rev(x,0)  =1/4
* pd_atan2_rev(x,-x) =3/8
* pd_atan2_rev(0,-x) =1/2
* pd_atan2_rev(-x,-x)=5/8
* pd_atan2_rev(-x,0) =3/4
* pd_atan2_rev(-x,x) =7/8

### sin of revolutions

	mlib: sin(x*(2.0*M_PI))
	C: pseudo_double pd_sin_rev(pseudo_double x);
	C++: PseudoDouble sin_rev(const PseudoDouble x);

### cos of revolutions

	mlib: cos(x*(2.0*M_PI))
	C: pseudo_double pd_cos_rev(pseudo_double x);
	C++: PseudoDouble cos_rev(const PseudoDouble x);

### atan2 in revolutions

	mlib: atan2(y,x)/(2.0*M_PI)
	C: pseudo_double pd_atan2_rev(pseudo_double y, pseudo_double x);
	C++: PseudoDouble atan2_rev(const PseudoDouble y, const PseudoDouble x);

### sin of radians

	mlib: sin(x)
	C: pseudo_double pd_sin(pseudo_double x);
	C++: PseudoDouble sin(const PseudoDouble x);

### cos of radians

	mlib: cos(x)
	C: pseudo_double pd_cos(pseudo_double x);
	C++: PseudoDouble cos(const PseudoDouble x);

### atan2 in radians

	mlib: atan2(y,x)
	C: pseudo_double pd_atan2(pseudo_double y, pseudo_double x);
	C++: PseudoDouble atan2(const PseudoDouble y, const PseudoDouble x);

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
