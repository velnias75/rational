/*
 * Copyright 2015 by Heiko Schäfer <heiko@rangun.de>
 *
 * This file is part of rational.
 *
 * rational is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * rational is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with rational.  If not, see <http://www.gnu.org/licenses/>.
 */

/**
 * @file
 *
 * This header contains specialization for use of
 * [the GNU Multiple Precision Arithmetic Library](https://gmplib.org/) as underlying
 * storage type.
 *
 * @author Heiko Schäfer <heiko@rangun.de>
 * @copyright 2015 by Heiko Schäfer <heiko@rangun.de>
 *
 * @defgroup gmp GNU Multiple Precisions extensions
 *
 * The header `gmp_rational.h` contains specializations especially for
 * [the GNU Multiple Precision Arithmetic Library](https://gmplib.org/) as underlying storage
 * type.\n
 * \n
 * If you use the *GMP extensions*, you'll need to link your application with `-lgmpxx -lgmp`
 */

#ifndef COMMONS_MATH_GMP_RATIONAL_H
#define COMMONS_MATH_GMP_RATIONAL_H

#include <gmpxx.h>

#include "rational.h"

#ifndef GMP_EPSILON
/**
 * @ingroup gmp
 * @def GMP_EPSILON
 *
 * The @c EPSILON used for approximating a float
 *
 * @see Commons::Math::EPSILON
 *
 * This define is passed @em as @em is to the contructor of @c mpz_class
 *
 * See [C++ Interface Integers]
 * (https://gmplib.org/manual/C_002b_002b-Interface-Integers.html#C_002b_002b-Interface-Integers)
 * for details
 */
#define GMP_EPSILON "1e-21", 30, 10
#endif

#if (__GNU_MP_VERSION * 10000 + __GNU_MP_VERSION_MINOR * 100 + __GNU_MP_VERSION_PATCHLEVEL) \
          < 50100
namespace std {

template<>
struct numeric_limits<mpz_class> {

    static const bool is_specialized = true;

    static mpz_class min() {
        return mpz_class();
    }

    static mpz_class max() {
        return mpz_class();
    }

    static mpz_class lowest() {
        return mpz_class();
    }

    static const int digits = 0;
    static const int digits10 = 0;
    static const int max_digits10 = 0;
    static const bool is_signed = true;
    static const bool is_integer = true;
    static const bool is_exact = true;
    static const int radix = 2;

    static mpz_class epsilon() {
        return mpz_class();
    }

    static mpz_class round_error() {
        return mpz_class();
    }

    static const int min_exponent = 0;
    static const int min_exponent10 = 0;
    static const int max_exponent = 0;
    static const int max_exponent10 = 0;
    static const bool has_infinity = false;
    static const bool has_quiet_NaN = false;
    static const bool has_signaling_NaN = false;
    static const float_denorm_style has_denorm = denorm_absent;
    static const bool has_denorm_loss = false;

    static mpz_class infinity() {
        return mpz_class();
    }

    static mpz_class quiet_NaN() {
        return mpz_class();
    }

    static mpz_class signaling_NaN() {
        return mpz_class();
    }

    static mpz_class denorm_min() {
        return mpz_class();
    }

    static const bool is_iec559 = false;
    static const bool is_bounded = false;
    static const bool is_modulo = false;
    static const bool traps = false;
    static const bool tinyness_before = false;
    static const float_round_style round_style = round_toward_zero;
};

}
#endif

namespace Commons {

namespace Math {

template<> struct TYPE_CONVERT<mpz_class> {

    inline explicit TYPE_CONVERT ( const mpz_class &v ) : val ( v ) {}

    template<typename U> inline U convert() const {
        return U ( val );
    }

private:
    const mpz_class &val;
};

template<> inline double TYPE_CONVERT<mpz_class>::convert<double>() const {
    return val.get_d();
}

template<> inline float TYPE_CONVERT<mpz_class>::convert<float>() const {
    return static_cast<float> ( val.get_d() );
}

template<> inline long signed int TYPE_CONVERT<mpz_class>::convert<long signed int>() const {
    return val.get_si();
}

template<> inline long unsigned int TYPE_CONVERT<mpz_class>::convert<long unsigned int>() const {
    return val.get_ui();
}
/**
 * @ingroup gmp
 * @ingroup gcd
 * @brief GMP GCD algorithm
 *
 * The gcd algorithm implemented in the GMP library
 *
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, bool IsSigned, template<class, typename = T, bool = IsSigned> class CHKOP,
         template<typename = T> class CONV> struct GCD_gmp;

template<template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_gmp<mpz_class, false, CHKOP, CONV> {

    inline mpz_class operator() ( const mpz_class &a, const mpz_class &b ) const {
        return GCD_gmp<mpz_class, true, CHKOP, CONV>() ( a, b );
    }

};

template<template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_gmp<mpz_class, true, CHKOP, CONV> {

    inline mpz_class operator() ( const mpz_class &a, const mpz_class &b ) const {

        mpz_class rop;

        mpz_gcd ( rop.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t() );

        return rop;
    }

};

template<> inline const mpf_class EPSILON<mpf_class>::value() {
    static mpf_class eps ( GMP_EPSILON );
    return eps;
}

template<>
inline mpf_class _approxFract<mpz_class, GCD_euclid, NO_OPERATOR_CHECK, mpf_class, true,
EPSILON, TYPE_CONVERT>::abs ( const mpf_class &nt ) const {
    return ::abs ( nt );
}

template<template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _lcm<mpz_class, GCD, CHKOP, true> {

    inline mpz_class operator() ( const mpz_class &a, const mpz_class &b ) const {

        mpz_class rop;

        mpz_lcm ( rop.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t() );

        return rop;
    }

};

/**
 * @ingroup gmp
 * @brief Rational class based on the GMP library
 */
typedef Rational<mpz_class, GCD_gmp, NO_OPERATOR_CHECK> gmp_rational;

}

}

#endif /* COMMONS_MATH_GMP_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
