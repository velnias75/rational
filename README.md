# rational
A C++ rational (fraction) template class

Include `rational.h` to be able to do fraction calculations. By simply including `rational.h` and 
specifying the storage type (any integer variant) you can create and use a fractional data type. 
For example, `Rational<long> foo(3, 4)` would create a fraction named `foo` with a value of `3/4`, 
and store the fraction using the `long` data type. 

Simple to use, fast, and modifiable for your needs, because its all in the one header file 
`rational.h` at your access. It includes all basic mathematical operations as well as comparison 
operators, and is very flexible. If you would like to see a feature implemented, just ask here.

The *storage type* should represent all integers within some (possibly infinite) interval 
including `0` and `1`. For example, the native `signed` or `unsigned int` and `long` types, or 
arbitrary-precision integers, may be used. Beyond ordinary integers, you should also be 
able to use any other [Euclidean domain](https://en.wikipedia.org/wiki/Euclidean_domain), 
perhaps not even an [ordered ring](https://en.wikipedia.org/wiki/Ordered_ring), but support 
for such types is experimental and has not been thoroughly tested. In fact, you should be able 
to use any [integral domain](https://en.wikipedia.org/wiki/Integral_domain), but you may need 
to apply a more sophisticated GCD algorithm; you can fall back to `GCD_null` if overflow is 
not a concern in practice. Finally, using non-integral domains is very likely to fail.

Features
--------

- exchangeable `GCD` algorithms *[`GCD_euclid_fast` (**default**), `GCD_euclid`and `GCD_stein`,
  choose as second template parameter, i. e. `Rational<long, GCD_stein> foo(3, 4)` for `Stein`]*
- optional *signed overflow/unsigned wrap* checking by throwing an `std::domain_error` exception
  (i.e. `Rational<storage_type, GCD_algo, Commons::Math::ENABLE_OVERFLOW_CHECK>`, default is
   no checking: `Commons::Math::NO_OVERFLOW_CHECK`)
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
- Support for 
    * [the GNU Multiple Precision Arithmetic Library](https://gmplib.org/) 
      (include `gmp_rational.h`)
    * [CLN - Class Library for Numbers](http://www.ginac.de/CLN/) (include `cln_rational.h`)

  as underlying storage type
- Expression templates for domain specific programming (include `expr_rational.h`)
- Construction of fractions from expression strings 
  (i.e. `Rational<long> expr("(11/2) * +(4.25+3.75)")`)
- Construction of fractions from continued fractions (from container of integer types)
- Extraction of continued fractions sequences from a fraction
- Construction of fractions from repeating decimals (i.e. `0.16666...` => `1/6`)

Notes for custom number types
-----------------------------

`rational.h` depends on following specialized fields of `std::numeric_limits<custom_type>`
- `std::numeric_limits<custom_type>::is_signed`
- `std::numeric_limits<custom_type>::is_integer`
- `std::numeric_limits<custom_type>::is_exact`

*(`gmp_rational.h` provides this specializations for GMP versions below 5.1)*

How to use
----------

See the test cases for examples on how to use `rational.h` with *C++ built-in types* as well
as with [the GNU Multiple Precision Arithmetic Library](https://gmplib.org/) or
the [CLN - Class Library for Numbers](http://www.ginac.de/CLN/) as underlying storage type.

* The header `gmp_rational.h` contains specializations especially for 
  [the GNU Multiple Precision Arithmetic Library](https://gmplib.org/) as underlying storage type.

  If you use the *GMP extensions*, you'll need to link your application with `-lgmpxx -lgmp`

* The header `cln_rational.h` contains specializations especially for the
  [CLN - Class Library for Numbers](http://www.ginac.de/CLN/) as underlying storage type.

  If you use the *CLN extensions*, you'll need to link your application with `-lcln`
