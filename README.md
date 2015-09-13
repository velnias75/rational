# rational
A C++ rational (fraction) template class

Include `rational.h` to be able to do fraction calculations. By simply including `rational.h` and 
specifying the storage type (any integer variant) you can create and use a fractional data type. 
For example, `Rational<long> foo(3, 4)` would create a fraction named `foo` with a value of `3/4`, 
and store the fraction using the `long` data type. 

Simple to use, fast, and modifiable for your needs, because its all in the one header file 
`rational.h` at your access. It includes all basic mathematical operations as well as comparison 
operators, and is very flexible. If you would like to see a feature implemented, just ask here.

**Whilst `Rational` tries to avoid overflows where possible, no overflow checking is performed
for performance reasons. Choose the storage type carefully and use `g++` features like `-ftrapv`
and `-fsanitize=undefined -fsanitize=signed-integer-overflow`.**

Features
--------

- exchangeable `GCD` algorithms *[`GCD_euclid` (**default**) and `GCD_stein`, choose as second
 template parameter, i. e. `Rational<long, GCD_stein> foo(3, 4)` for `Stein`]*
- optimized for `signed` and `unsigned` types
- additional operators: 
  - `mod` to split inproper fractions in integer and fraction part
  - `abs` to get the absolute value
  - `%` (modulo) for rational modulo
  - increment (`++x` & `x++`) and decrement (`--x` & `x--`)
  - unary `plus`and `minus`
- Construction of inproper (mixed) fractions, i.e. `Rational<long> foo(2, 3, 4)`for `2 3/4` resp. 
  `2.75`
- Construction of approximate fractions, i.e. `Rational<long> foo(3.14159265358979323846)`for `π` 
  resp. `245850922/78256779` *(approximation is dependent on compiler and chosen storage type)*

Notes for custom number types
-----------------------------

`rational.h` depends on following specialized fields of `std::numeric_limits<custom_type>`
- `std::numeric_limits<custom_type>::is_signed`
- `std::numeric_limits<custom_type>::is_integer`
- `std::numeric_limits<custom_type>::is_exact`
