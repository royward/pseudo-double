# pseudo-double

A relatively fast C and C++ 64 bit floating point library written using only integer operations for cross platform consistency.

# Overview

This is a partial floating point library that only uses integer CPU instructions (no floating point) and is designed to provide consistent cross platform results in cases where fixed point does not have enough dynamic range. Exact consistency across platforms is particularly important for procedural generation or other very ill conditioned systems where a tiny change in the input or calculation can produce an enormous change in the outputs. It was inspired by this post:

https://gamedev.stackexchange.com/questions/202475/consistent-cross-platform-procedural-generation

There are many libraries that exist to provide floating point functionality on integer processors. They all (that I have found) follow the IEEE 754 standard at least in part, and the vast majority of them are only single precision. The pseudo-double library does not follow IEEE 754 standard at all, and instead has design choices more suited to a software implementation. This results in a library that (on x86-64) runs 4-15 times slower than the hardware floating point.

This library has both C and C++ bindings, and it has been tested on x86-x64 with gcc/g++ and ARMv8-A with clang.

It also allows the number of expononent bits to be set with a macro, which enables choices to be made in the range/precision tradeoff.

# Use cases

The main use cases for this code are:

* Something where cross-platform consistency is important and there is too much dynamic range required to use fixed point.

* Applications where a range/precision tradeoff different than the IEEE 754 one is useful setting (see PSEUDO_DOUBLE_EXP_BITS below).

* This may be adaptable to hardware that doesn't support floating point but double precision or the speedup of not requiring IEEE 754 is still desired. The would likely require considerable rework, as a processor that doesn't already support floating point is unlikely to support things like 128 bit integers.

# Usage

See Functions.md for details of all the provided functions.

For most uses, the easiest option is to include the source files directly in your project.

## C

pseudo_double is a struct containing a single number to provide type safety.

If you want to avoid use of the struct, pseudo_double_i aliased to uint64_t, and there there are functions pdi\_\*() and \*\_pdi() that match the pd\_\*() and \*\_pd() functions. Use at your own risk.

### Files

**pseudo_double.h**: file to include to access the C functionality. Some of the basic functions are provided by inline functions.

**pseudo_double.c**: Add this into your project. This includes the functionality not provided by pseudo_double.h

## C++

PseudoDouble is provided as a class with operator overloading and function overloading to provide the same syntax as the standard floating point operations.

### Files

**PseudoDouble.h**: file to include to access the C++ functionality

**pseudo_double.cpp**: is just a C++ wrapper for pseudo_double.c, or pseudo_double.c can be used directly.

**pseudo_double.h**: is included by PseudoDouble.h.

**PseudoDouble_test.cpp**: a test for the rest of the library.

**PseudoDouble_speed_test.cpp**: a simple speed test. 10x10 matrix inversion plus some short loop tests.

## Example code

Find the roots of $0.3x^2-4x+6$ using the quadratic formula $\frac{-b\pm\sqrt{b^2-4ac}}{2a}$

### C/C++ With doubles

This is some code with doubles as it might normally be done.

	double a=0.3;
	double b=-4.0;
	double c=6.0;
	double disc=sqrt(b*b-4.0*a*c);
	double sol1=(-b-disc)/(2.0*a);
	double sol2=(-b+disc)/(2.0*a);
	printf("C: Solution 1 = %lf\n",sol1);
	printf("C: Solution 2 = %lf\n",sol2);

### C with pseudo_double

	pseudo_double a=int64fixed10_to_pd(3,-1); // 0.3
	pseudo_double b=int64_to_pd(-4);
	pseudo_double c=int64_to_pd(6);
	pseudo_double disc=pd_sqrt(pd_sub(pd_mult(b,b),pd_mult(pd_mult(int64_to_pd(4),a),c)));
	pseudo_double sol1=pd_div(pd_sub(pd_neg(b),disc),pd_mult(int64_to_pd(2),a));
	pseudo_double sol2=pd_div(pd_add(pd_neg(b),disc),pd_mult(int64_to_pd(2),a));
	printf("C: Solution 1 = %lf\n",pd_to_double(sol1));
	printf("C: Solution 2 = %lf\n",pd_to_double(sol2));

Note that we have to be careful constructing the a=0.3 The code:

	pseudo_double a=double_to_pd(0.3);

would also work, but that is dependent on hardware floating point and may not be guaranteed to give the same result cross platform.

### C++ with PseudoDouble

Here we get to make good use of operator and function overloading.

	PseudoDouble a=PD_create_fixed10(3,-1); // 0.3
	PseudoDouble b=-4;
	PseudoDouble c=6;
	PseudoDouble disc=sqrt(b*b-PseudoDouble(4)*a*c);
	PseudoDouble sol1=(-b-disc)/(PseudoDouble(2)*a);
	PseudoDouble sol2=(-b+disc)/(PseudoDouble(2)*a);
	std::cout << "C++: Solution 1 = " << sol1 << std::endl;
	std::cout << "C++: Solution 2 = " << sol2 << std::endl;

### C with pseudo_double_i (type unsafe)

Here we don't use the pseudo_double struct but use the pseudo_double_i instead. It is very easy to accidentally use direct integer operations and get garbage. Not recommended.

	pseudo_double_i a=int64fixed10_to_pdi(3,-1); // 0.3
	pseudo_double_i b=int64_to_pdi(-4);
	pseudo_double_i c=int64_to_pdi(6);
	pseudo_double_i disc=pdi_sqrt(pdi_sub(pdi_mult(b,b),pdi_mult(pdi_mult(int64_to_pdi(4),a),c)));
	pseudo_double_i sol1=pdi_div(pdi_sub(pdi_neg(b),disc),pdi_mult(int64_to_pdi(2),a));
	pseudo_double_i sol2=pdi_div(pdi_add(pdi_neg(b),disc),pdi_mult(int64_to_pdi(2),a));
	printf("C (unsafe): Solution 1 = %lf\n",pdi_to_double(sol1));
	printf("C (unsafe): Solution 2 = %lf\n",pdi_to_double(sol2));

## Running the tests

To built the library test:

	g++ -Wall -Wextra ${ANY_OTHER_FLAGS} -o pseudo_double_test pseudo_double.cpp PseudoDouble_test.cpp

When the generated executable is run, it will print out a line for any failed tests, followed by a count of tests passed/done:

	./pseudo_double_test 
	Tests done, passed 3832004/3832004

Similarly, the speed tests can be built using:

	g++ -Wall -Wextra ${ANY_OTHER_FLAGS} -o pseudo_double_speed_test pseudo_double.cpp PseudoDouble_speed_test.cpp

# Funtions supported

See Functions.md for details of all the provided functions.

* **Basic operators:** +, - (both unary and binary), *, /, ==, !=, >, >=, <, <=

* **Functions matching <math.h>**: floor, ceil, round, sqrt, ldexp, exp2, exp, log2, log, log10, pow, sin, cos, atan2, abs, fabs

* **Functions not found in <math.h>**: inv_sqrt, sin_rev, cos_rev, atan2_rev, conversion to and from doubles, pseudo-double creation

* **Functions not currently supported by pseudo-double**:  acos, asin, tan, atan, hyperbolic trigonometry, frexp, expm1, ilogb, log1p, logb, scalbn, scalbln, cbrt, hypot, erf, erfc, tgamma, lgamma, fmod, trunc, lround, llround, rint, lrint, llrint, nearbyint, remainder, remquot, copysign, nan, nextafter, nexttoward, fdim, fmax, fmin, fma, fpclassify, signbit, isfinite, isinf, isnan, isnormal, all the comparison macros.

# Overflows

There is an error value (PD_NAN in C) overflows/errors/out_of_range are represented by all bits set to 1.

There are three macros that can be set to determine the behaviour in these cases:

PD_DO_ERROR_OVERFLOW  
C default: return PD_NAN  
C++ default: throw std::overflow_error("overflow")

PD_DO_ERROR_UNDERFLOW  
C default: return 0  
C++ default: return PseudoDouble(0)

PD_DO_ERROR_RANGE  
C default: return PD_NAN  
C++ default: throw std::range_error("range")

NOTE: in the interests of performance, PD_NAN is _not_ checked for on input, so if it is used, it should be checked for explicitly after any calculation that might generate that value.

# Other settings

## PSEUDO_DOUBLE_EXP_BITS

This sets the number of bits in the exponent, defaulting to 16 if not set.

If the number of the exponent bits is $n$:

* The mantissa will have $64-n$ fits for a precision of $63-n$.

* Any number $x$ will be in the range $-2^{2^{n-1}-2} \le x < -2^{-2^{n-1}-2}$ or $2^{-2^{n-1}-2} \le x < 2^{2^{n-1}-2}$

Setting PSEUDO_DOUBLE_EXP_BITS to 8, 16 or 32 will be slightly faster than other values, as the CPU can take advantage of the internal integer sizes rather than having to do shifts.

## PD_ERROR_CHECK

Overflow, range and some underflow checking can be turned off by setting the macro PD_ERROR_CHECK to 0 (default is 1). This may give a very slight preformance increase, but at the cost of returning undetectable garbage instread of and error. It is not worth turning errors off unless you are certain that overflow/range/underflow errors will not occur. This will also cause some "may be used uninitialized in this function" errors on compilation.

# Design Considerations

This library is designed to be a tradeoff between speed and accuracy. It does not get full IEEE 754 accuracy although it is often close, but should be reasonably performant, although of course not even close to native floating point.

Four properties to consider when determining how to perform calculations on continuous quantities (things that would be represented mathematically with real numbers): precision, speed, (dynamic) range and (cross platform) consistency.

| Type          | Size (bits) | Exponent (bits) | Precision (bits) | Range                                                    | Time | Consistency |
| ------------- | ----------- | --------------- | ---------------- | -------------------------------------------------------- | -----| ----------- |
| float         | 32          | 8               | 53               | $\pm(1.18\times 10^{-38} \ldots 3.4\times 10^{38})$      | 1    | no          |
| double        | 64          | 11              | 24               | $\pm(2.23\times 10^{-308} \ldots 1.8\times 10^{308})$    | 1    | no          |
| fixed         | 32          | 0               | 0 ... 31 $\ast$  | $\pm c(1 \ldots 2.1\times 10^{9})$                       | 1    | yes         |
| fixed         | 64          | 0               | 0 ... 63 $\ast$  | $\pm c(1 \ldots 9.2\times 10^{18})$                      | 1    | yes         |
| pseudo-double | 64          | 16 $\dagger$    | 47               | $\pm(2.10\times 10^{-4933} \ldots 2.97\times 10^{4931})$ | ~10  | yes         |
| pseudo-double | 64          | 8 $\dagger$     | 55               | $\pm(7.35\times 10^{-40} \ldots 8.51\times 10^{37})$     | ~10  | yes         |

$\ast$ precision for fixed depends on how much of the 32-bit or 64 bit word is used.

$c$ for fixed point numbers is a positive constant, typically of the form $2^{-k}, k \ge 0$

$\dagger$ these are representative values, exponent can be any value from 1..62, with corresponding changes in precision and range. 

Pseudo Doubles are intended to give precision, range and consistency while sacrificing as little speed as possible, although they will never be as fast as float, double or fixed.

# Underlying structure

Because this runs on a CPU using integer instructions rather than dedicated hardware, different choices are made for the storage format than IEE 754.

For reference, an IEEE 754 double has a sign bit, 11 bits of exponent, and a 52 bit unsigned mantissa with the most significant bit removed:

	seeeeeeeeeeemmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
	6666555555555544444444443333333333222222222211111111110000000000
	3210987654321098765432109876543210987654321098765432109876543210

0 is represented by all zeros, -0 is represented by a 1 sign bit, the rest zero.

Pseudo doubles have a 48 bit signed mantissa (no bits removed), and a settable (default 16 bit) exponent:

	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmeeeeeeeeeeeeeeee
	6666555555555544444444443333333333222222222211111111110000000000
	3210987654321098765432109876543210987654321098765432109876543210

Here is the layout with an 8 bit exponent:

	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmeeeeeeee
	6666555555555544444444443333333333222222222211111111110000000000
	3210987654321098765432109876543210987654321098765432109876543210


A number with b exponent bits is represented as $m 2^{e-2^{b-1}}$ where m is a signed fixed point, $-0.5 \le m < 0.5$ . For instance, with a 16 bit exponent, a number is represented by $m 2^{e-32768}$ .

0 is represented by all zeros, there is no -0.

There is an overflow value represented by all 1 bits. It is is returned in cases where an overflow, infinity or error happens.

* 48 (or 56) bits of mantissa rather than 52 bits means everything is on 8 bit boundaries, which is a slight performance improvement.
* a signed mantissa rather than a separate sign bit means less branching and branches are more predictable.
* putting the mantissa in the most significant bits means that comparison with less than, greater than, or equal to zero can be does by looking at the integer values.
* the exponent in the least significant bits means that ldexp can be done with increment/decrement (except in the case of overflow).
* pseudo doubles don't use denormalized representations - it's a lot of extra checking (cheap in hardware, expensive in software) for only a small increase in dynamic range. Unless set up to return some sort of error, underflows go straight to zero.

Except for zero and PD_NAN, the most significant bits of the mantissa after normalization are always different (01 or 10). Technically this redundancy could be removed to give one extra bit of precision, but that would result in more complex computation.

There is an asymmetry that needs to be considered caused by $-0.5 \le m < 5$. Most numbers can be negated by just leaving the exponent the same and integer negating the mantissa. The exception is when the number represents a power of two or negative power of two:

A power of two is represented by:

	01000000.... * 2^exp

Straight negating this would give:

	11000000.... * 2^exp

This then need to be normalized to give:

	10000000.... * 2^(exp-1)

This process is reversed for negating a negative power of two.

The upside is that because this case only occurs with power of two or negative power of two, this case is likely to be rarer than the other cases, which will improve branch prediction.

# Porting to other compilers

There are two non-standard features that are gcc/g++/clang specific and might need to be adjusted for other compilers:

### 128 bit signed integers

gcc/g++/clang uses __int128. This is required for multiplication, division, pow and atan2.

### Count leading zeros

gcc/g++/clang uses __builtin_clzll. This is required for all functions. For x86-64 the intel intrinsic _lzcnt_u64 is available. 

# Other notes

The constants used for generating transcendental functions were generated using lolremez (https://github.com/samhocevar/lolremez) which is an implementation of the Remez algorithm (https://en.wikipedia.org/wiki/Remez_algorithm). the lolremez commands and results are included as comments in the code. Adjustments were then made:
* scaled integers are used rather than floating point
* in some cases the calculation has various scale changes applied to best use as much of the range of 64 bit integers as possible without overflowing
* small adjustments are made in the last step to try and ensure the overall function is as smooth as possible at the boundaries, possibly at the expense of some accuracy.
