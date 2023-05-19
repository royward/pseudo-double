# pseudo-float

A relatively fast C and C++ 64 bit floating point library written using only integer operations for cross platform consistency.

# Overview

This is a partial floating point library that only uses integer CPU instructions (no floating point) and is designed to provide consistent cross platform results in cases where fixed point does not have enough dynamic range. Exact consistency across platforms is particularly important for procedural generation or other very ill conditioned systems where a tiny change in the input or calculation can produce an enormous change in the outputs. It was inspired by this post:

https://gamedev.stackexchange.com/questions/202475/consistent-cross-platform-procedural-generation

There are many libraries that exist to provide floating point functionality on integer processors. They all (that I have found) follow the IEEE 754 standard at least in part, and the vast majority of them are only single precision. The pseudo-float library does not follow IEEE 754 standard at all, and instead has design choices more suited to a software implementation. This results in a library that (on x86-64) runs 4-15 times slower than the hardware floating point.

This library has both C and C++ bindings, and it has been tested on x86-x64 with gcc/g++ and ARMv8-A with clang. 

# Usage

See Functions.md for details of all the provided functions.

For most uses, the easiest option is to include the source files directly in your project.

## C

pseudo_float is aliased to uint64_t

### Files

**pseudo_float.h**: file to include to access the C functionality. Some of the basic functions are provided by inline functions.

**pseudo_float.c**: Add this into your project. This includes the functionality not provided by pseudo_float.h

## C++

PseudoFloat is provided as a class with operator overloading and function overloading are used to provide the same syntax as the standard floating point operations.

### Files

**PseudoFloat.h**: file to include to access the C++ functionality

**pseudo_float.cpp**: is just a C++ wrapper for pseudo_float.c

**PseudoFloat_test.cpp**: a test for the rest of the library.

**PseudoFloat_speed_test.cpp**: a simple speed test. 10x10 matrix inversion plus some short loop tests.

### Running the tests

To built the library test:

	g++ -Wall -Wextra ${ANY_OTHER_FLAGS} -o pseudo_float_test pseudo_float.cpp PseudoFloat_test.cpp

When the generated executable is run, it will print out a line for any failed tests, followed by a count of tests passed/done:

	./pseudo_float_test 
	Tests done, passed 3863863/3863863

Similarly, the speed tests can be built using:

	g++ -Wall -Wextra ${ANY_OTHER_FLAGS} -o pseudo_float_speed_test pseudo_float.cpp PseudoFloat_speed_test.cpp

# Funtions supported

See Functions.md for details of all the provided functions.

* **Basic operators:** +, - (both unary and binary), *, /, ==, !=, >, >=, <, <=

* **Functions matching <math.h>**: floor, ceil, round, sqrt, ldexp, exp2, exp, log2, log, log10, pow, sin, cos, atan2, abs, fabs

* **Functions not found in <math.h>**: inv_sqrt, sin_rev, cos_rev, atan2_rev, conversion to and from doubles, pseudo-float creation

* **Functions not currently supported by pseudo-float**:  acos, asin, tan, atan, hyperbolic trigonometry, frexp, expm1, ilogb, log1p, logb, scalbn, scalbln, cbrt, hypot, erf, erfc, tgamma, lgamma, fmod, trunc, lround, llround, rint, lrint, llrint, nearbyint, remainder, remquot, copysign, nan, nextafter, nexttoward, fdim, fmax, fmin, fma, fpclassify, signbit, isfinite, isinf, isnan, isnormal, all the comparison macros.

# Overflows

There is an error value (PF_NAN in C) overflows/errors/out_of_range are represented by all bits set to 1.

There are three macros that can be set to determine the behaviour in these cases:

PF_DO_ERROR_OVERFLOW  
C default: return PF_NAN  
C++ default: throw std::overflow_error("overflow")

PF_DO_ERROR_UNDERFLOW  
C default: return 0  
C++ default: return PseudoFloat(0)

PF_DO_ERROR_RANGE  
C default: return PF_NAN  
C++ default: throw std::range_error("range")

NOTE: in the interests of performance, PF_NAN is _not_ checked for on input, so if it is used, it should be checked for explicitly after any calculation that might generate that value.

Overflow, range and some underflow checking can be turned off by setting the macro PF_ERROR_CHECK to 0 (default is 1). This may give a very slight preformance increase, but at the cost of returning undetectable garbage instread of and error. It is not worth turning errors off unless you are certain that overflow/range/underflow errors will not occur. This will also cause some "may be used uninitialized in this function" errors on compilation.

# Design Considerations

This library is designed to be a tradeoff between speed and accuracy. It does not get full IEEE 754 accuracy although it is often close, but should be reasonably performant, although of course not even close to native floating point.

Four properties to consider when determining how to perform calculations on continuous quantities (things that would be represented mathematically with real numbers): precision, speed, (dynamic) range and (cross platform) consistency.

	float:  32             speed  range
	double: 64  precision  speed  range
	fixed:  32             speed         consistency
	fixed:  64  precision  speed         consistency
	pseudo: 64  precision         range  consistency

Pseudo Floats are intended to give precision, range and consistency while sacrificing as little speed as possible, although they will never be as fast as float, double or fixed. The pseudo-float library does not use 

# Underlying structure

Because this runs on a CPU using integer instructions rather than dedicated hardware, different choices are made for the storage format than IEE 754.

For reference, an IEEE 754 double has a sign bit, 11 bits of exponent, and a 52 bit unsigned mantissa with the most significant bit removed:

	seeeeeeeeeeemmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmm
	6666555555555544444444443333333333222222222211111111110000000000
	3210987654321098765432109876543210987654321098765432109876543210

0 is represented by all zeros, -0 is represented by a 1 sign bit, the rest zero.

Pseudo floats have a 48 bit signed mantissa (no bits removed), and a 16 bit exponent:

	mmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmmeeeeeeeeeeeeeeee
	6666555555555544444444443333333333222222222211111111110000000000
	3210987654321098765432109876543210987654321098765432109876543210

A number is represented as m*2^(e-16384) where m is a signed fixed point, -0.5<=m<0.5
0 is represented by all zeros, there is no -0.
There is an overflow value represented by all 1 bits. It is is returned in cases where an overflow, infinity or error happens.

* 48 bits of mantissa rather than 52 bits means everything is on 8 bit boundaries, which is a slight performance improvement.
* a signed mantissa rather than a separate sign bit means less branching and branches are more predictable.
* putting the mantissa in the most significant bits means that comparison with less than, greater than, or equal to zero can be does by looking at the integer values.
* the exponent in the least significant bits means that ldexp can be done with increment/decrement (except in the case of overflow).
* pseudo floats don't use denormalized representations - it's a lot of extra checking (cheap in hardware, expensive in software) for only a small increase in dynamic range. Unless set up to return some sort of error, underflows go straight to zero.

Except for zero and PF_NAN, the most significant bits of the mantissa after normalization are always different (01 or 10). Technically this redundancy could be removed to give one extra bit of precision, but that would result in more complex computation.

There is an asymmetry that needs to be considered caused by -0.5<=m<5. Most numbers can be negated by just leaving the exponent the same and integer negating the mantissa. The exception is when the number represents a power of two or negative power of two:

A power of two is represented by:

	01000000.... * 2^exp

Straight negating this would give:

	11000000.... * 2^exp

This then need to be normalized to give:

	10000000.... * 2^(exp-1)

This process is reversed for negating a negative power of two.

The upside is that because this case only occurs with power of two or negative power of two, this case is likely to be rarer than the other cases, which will improve branch prediction.

# Other notes

The constants used for generating transcendental functions were generated using lolremez (https://github.com/samhocevar/lolremez) which is an implementation of the Remez algorithm (https://en.wikipedia.org/wiki/Remez_algorithm). the lolremez commands and results are included as comments in the code. Adjustments were then made:
* scaled integers are used rather than floating point
* in some cases the calculation has various scale changes applied to best use as much of the range of 64 bit integers as possible without overflowing
* small adjustments are made in the last step to try and ensure the overall function is as smooth as possible at the boundaries, possibly at the expense of some accuracy.


