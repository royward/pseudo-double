# Change Log

# 1.1.0 - 2024-03-03

* Fixes to ldexp(0,n)
* Fix for int64fixed10_to_pdi / int64fixed10_to_pd / PD_create_fixed10 for negative numbers
* NEW: Rust binding

# 1.0.6 - 2024-03-24

* Added a namespace 'pseudodouble' for the C++ version
* Added functions for max and min
* fixed a missing <string> include for the C++ version
* Small documentation fixes

# Change Log

# 1.0.5 - 2023-12-10

* Fixed bug where uint64_to_pdi treated argument as signed rather than unsigned

# 1.0.4 - 2023-11-19

* Fixed the zero comparison operators 

# 1.0.3 - 2023-08-22

* Added support for Visual C++
* Added C compare with zero functions
  
# 1.0.2 - 2023-05-24

* Fixed some C++isms that crept into the C code

# 1.0.1 - 2023-05-20

* Changed **pseudo_double** from a typedef to a struct to assist with type safety
* **pseudo_double_i** is now the typedef

# 1.0.0 - 2023-05-20

* Changed name from pseudo-float to pseudo-double
* First public release
