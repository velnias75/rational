/*
 * Copyright 2015-2016 by Heiko Schäfer <heiko@rangun.de>
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
 * @copyright 2015-2016 by Heiko Schäfer <heiko@rangun.de>
 *
 * @defgroup gmp GNU Multiple Precisions extensions
 *
 * The header `gmp_rational.h` contains specializations especially for
 * [the GNU Multiple Precision Arithmetic Library](https://gmplib.org/)
 * as underlying storage type.\n
 * \n If you use the *GMP extensions*, you'll need to link your
 * application with `-lgmpxx -lgmp`
 */

#ifndef COMMONS_MATH_GMP_RATIONAL_H
#define COMMONS_MATH_GMP_RATIONAL_H

#if defined(IN_IDE_PARSER)
#define HAVE_MPREAL_H 1
#endif

#if (_POSIX_C_SOURCE >= 200809L || _XOPEN_SOURCE >= 700) || _GNU_SOURCE
#include <cstring>
#endif

#include <gmpxx.h>

#ifdef HAVE_MPREAL_H
#include <mpreal.h>
#endif

#include "rational.h"

#ifndef GMP_EPSILON
/**
 * @ingroup gmp
 * @def GMP_EPSILON
 *
 * @brief The @c EPSILON used for approximating a float
 *
 * @see Commons::Math::EPSILON
 *
 * This define is passed @em as @em is to the contructor of @c mpz_class
 *
 * See [C++ Interface Integers]
 * (https://gmplib.org/manual/C_002b_002b-Interface-Integers.html#C_002b_002b-Interface-Integers)
 * for details
 */
#define GMP_EPSILON "1e-100"
#endif

namespace std {

template<> inline void swap<mpz_class> ( mpz_class &x, mpz_class &y ) {
    mpz_swap ( x.get_mpz_t(), y.get_mpz_t() );
}

#if (defined(__GMP_MP_RELEASE) && __GMP_MP_RELEASE < 50103) || \
     (defined(__GNU_MP_RELEASE) && __GNU_MP_RELEASE < 50103)
template<> struct numeric_limits<mpz_class> {

    static const bool is_specialized = true;

    static const mpz_class min() {
        return mpz_class();
    }

    static const mpz_class max() {
        return mpz_class();
    }

    static const mpz_class lowest() {
        return mpz_class();
    }

    static const int digits = 0;
    static const int digits10 = 0;
    static const int max_digits10 = 0;
    static const bool is_signed = true;
    static const bool is_integer = true;
    static const bool is_exact = true;
    static const int radix = 2;

    static const mpz_class epsilon() {
        return mpz_class();
    }

    static const mpz_class round_error() {
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

    static const mpz_class infinity() {
        return mpz_class();
    }

    static const mpz_class quiet_NaN() {
        return mpz_class();
    }

    static const mpz_class signaling_NaN() {
        return mpz_class();
    }

    static const mpz_class denorm_min() {
        return mpz_class();
    }

    static const bool is_iec559 = false;
    static const bool is_bounded = false;
    static const bool is_modulo = false;
    static const bool traps = false;
    static const bool tinyness_before = false;
    static const float_round_style round_style = round_toward_zero;
};
#endif

#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
template<> struct not_equal_to<mpz_class> : public binary_function<mpz_class, mpz_class, bool> {
    bool operator() ( const mpz_class &x, const mpz_class &y ) const {
        return mpz_cmp ( x.get_mpz_t(), y.get_mpz_t() );
    }
};
#pragma GCC diagnostic pop

}

namespace Commons {

namespace Math {

template<> struct ExpressionEvalTraits<mpz_class> {
#ifdef HAVE_MPREAL_H
    typedef mpfr::mpreal NumberType;
#else
    typedef mpf_class NumberType;
#endif
};

template<> struct DecomposeBaseTraits<mpz_class, true> {
    typedef mpz_class digit_type;
    enum { Base = 10 };
};

template<> struct _type_round_helper<mpz_class> {
    mpz_class operator() ( const mpz_class &tr ) const {
        return tr;
    }
};

template<> inline mpz_class TYPE_CONVERT<long double>::convert<mpz_class>() const {
    std::ostringstream os;
    os.precision ( std::numeric_limits<double>::digits );
    os << static_cast<double> ( val );
    // Note: gmp can handle only doubles, so we need to narrow it to double
    return mpz_class ( os.str() );
}

template<> inline mpf_class TYPE_CONVERT<std::string>::convert<mpf_class>() const {
    return mpf_class ( val.c_str() );
}

#if (_POSIX_C_SOURCE >= 200809L || _XOPEN_SOURCE >= 700) || _GNU_SOURCE
class _gmp_charptr_convert_helper {

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
public:
    _gmp_charptr_convert_helper ( const _gmp_charptr_convert_helper & ) = delete;
    _gmp_charptr_convert_helper &operator= ( const _gmp_charptr_convert_helper & ) = delete;
#else
    _gmp_charptr_convert_helper ( const _gmp_charptr_convert_helper & );
    _gmp_charptr_convert_helper &operator= ( const _gmp_charptr_convert_helper & );

public:
#endif

    _gmp_charptr_convert_helper ( const char *val, const char *len )
        : val_ ( strndup ( val, static_cast<std::size_t> ( len - val ) ) ) {}

    ~_gmp_charptr_convert_helper() {
        std::free ( val_ );
    }

    char *c_str() const {
        return val_;
    }

private:
    char *val_;
};
#endif

template<> inline mpf_class TYPE_CONVERT<const char *>::convert<mpf_class>() const {
#if (_POSIX_C_SOURCE >= 200809L || _XOPEN_SOURCE >= 700) || _GNU_SOURCE
    return mpf_class ( len ? _gmp_charptr_convert_helper ( val, len ).c_str() : val );
#else
    return mpf_class ( len ? std::string ( val, len ).c_str() : val );
#endif
}

#ifdef HAVE_MPREAL_H
template<> inline mpfr::mpreal TYPE_CONVERT<std::string>::convert<mpfr::mpreal>() const {
    return mpfr::mpreal ( val.c_str() );
}

template<> inline mpfr::mpreal TYPE_CONVERT<const char *>::convert<mpfr::mpreal>() const {
    return mpfr::mpreal ( len ? std::string ( val, len ) : val );
}
#endif

template<> struct TYPE_CONVERT<mpz_class> {

    RATIONAL_NOCOPYASSIGN ( TYPE_CONVERT<mpz_class> );

    explicit TYPE_CONVERT ( const mpz_class &v ) : val ( v ) {}

    template<typename U> U convert() const;

private:
    const mpz_class &val;
};

template<typename U> U TYPE_CONVERT<mpz_class>::convert() const {
    return U ( val );
}

template<> inline double TYPE_CONVERT<mpz_class>::convert<double>() const {
    return val.get_d();
}

template<> inline long double TYPE_CONVERT<mpz_class>::convert<long double>() const {
    return static_cast<long double> ( val.get_d() );
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

#ifdef HAVE_MPREAL_H
template<> inline mpfr::mpreal TYPE_CONVERT<mpz_class>::convert<mpfr::mpreal>() const {
    return val.get_mpz_t();
}

template<> struct TYPE_CONVERT<mpfr::mpreal> {

    RATIONAL_NOCOPYASSIGN ( TYPE_CONVERT<mpfr::mpreal> );

    explicit TYPE_CONVERT ( const mpfr::mpreal &v ) : val ( v ) {}

    template<typename U> U convert() const {

        U out;

        mpfr_get_z ( out.get_mpz_t(), val.mpfr_ptr(), val.get_default_rnd() );

        return out;
    }

private:
    const mpfr::mpreal &val;
};
#endif

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

    mpz_class operator() ( const mpz_class &a, const mpz_class &b ) const {
        return GCD_gmp<mpz_class, true, CHKOP, CONV>() ( a, b );
    }

};

template<template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_gmp<mpz_class, true, CHKOP, CONV> {

    mpz_class operator() ( const mpz_class &a, const mpz_class &b ) const {
#if (defined(__GMP_MP_RELEASE) && __GMP_MP_RELEASE < 60100) || \
     (defined(__GNU_MP_RELEASE) && __GNU_MP_RELEASE < 60100)
        mpz_class rop;

        mpz_gcd ( rop.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t() );

        return rop;
#else
        return gcd ( a, b );
#endif
    }

};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc> struct _abs<mpz_class, GCD, CHKOP, Alloc, true> {

    Rational<mpz_class, GCD, CHKOP, Alloc> operator()
    ( const Rational<mpz_class, GCD, CHKOP, Alloc> &r ) const {
        return Rational<mpz_class, GCD, CHKOP, Alloc> ( abs ( r.numerator() ), r.denominator() );
    }
};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc> struct _abs<mpz_class, GCD, CHKOP, Alloc, false> {

    Rational<mpz_class, GCD, CHKOP, Alloc> operator()
    ( const Rational<mpz_class, GCD, CHKOP, Alloc> &r ) const {
        return _abs<mpz_class, GCD, CHKOP, Alloc, true>() ( r );
    }
};

template<> inline const mpf_class EPSILON<mpf_class>::value() {
    static mpf_class eps ( GMP_EPSILON );
    return eps;
}

#ifdef HAVE_MPREAL_H
template<> inline const mpfr::mpreal EPSILON<mpfr::mpreal>::value() {
    static mpfr::mpreal eps ( GMP_EPSILON );
    return eps;
}

template<template<typename> class EPSILON> struct _approxUtils<mpfr::mpreal, EPSILON> {

    static bool approximated ( const mpfr::mpreal &af, const mpfr::mpreal &nt ) {
        return !mpfr_cmp ( af.mpfr_ptr(), nt.mpfr_ptr() );
    }

    static mpfr::mpreal reciprocal ( const mpfr::mpreal &x ) {
        return 1u / x;
    }

    const static mpfr::mpreal eps_;

private:
    static mpfr::mpreal abs ( const mpfr::mpreal &nt ) {
        return mpfr::abs ( nt );
    }
};

template<template<typename> class EPSILON>
const mpfr::mpreal _approxUtils<mpfr::mpreal, EPSILON>::eps_ ( EPSILON<mpfr::mpreal>::value() );
#endif

template<template<typename> class EPSILON> struct _approxUtils<mpf_class, EPSILON> {

    static bool approximated ( const mpf_class &af, const mpf_class &nt ) {
        return mpf_eq ( af.get_mpf_t(), nt.get_mpf_t(), std::min ( af.get_prec(), nt.get_prec() ) );
    }

    static mpf_class reciprocal ( const mpf_class &x ) {
        mpf_class rop ( 1.0, x.get_prec() );
        mpf_div ( rop.get_mpf_t(), rop.get_mpf_t(), x.get_mpf_t() );
        return rop;
    }

    const static mpf_class eps_;

private:
    static mpf_class abs ( const mpf_class &nt ) {
        return ::abs ( nt );
    }
};

template<template<typename> class EPSILON>
const mpf_class _approxUtils<mpf_class, EPSILON>::eps_ ( EPSILON<mpf_class>::value() );

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc> struct _swapSign<mpz_class, GCD, CHKOP, Alloc, true> {

    Rational<mpz_class, GCD, CHKOP, Alloc> &
    operator() ( Rational<mpz_class, GCD, CHKOP, Alloc> &r ) const {

        if ( sgn ( r.m_denom ) == -1 ) {
            r.m_numer = -r.m_numer;
            r.m_denom = -r.m_denom;
        }

        return r;
    }

};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc> struct _lcm<mpz_class, GCD, CHKOP, Alloc, true> {

    mpz_class operator() ( const mpz_class &a, const mpz_class &b ) const {
#if (defined(__GMP_MP_RELEASE) && __GMP_MP_RELEASE < 60100) || \
     (defined(__GNU_MP_RELEASE) && __GNU_MP_RELEASE < 60100)
        mpz_class rop;

        mpz_lcm ( rop.get_mpz_t(), a.get_mpz_t(), b.get_mpz_t() );

        return rop;
#else
        return lcm ( a, b );
#endif
    }

};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc> struct _psq<mpz_class, GCD, CHKOP, Alloc> {

    Rational<mpz_class, GCD, CHKOP, Alloc> operator()
    ( const Rational<mpz_class, GCD, CHKOP, Alloc> &x,
      const Rational<mpz_class, GCD, CHKOP, Alloc> &y ) const {

        const typename Rational<mpz_class, GCD, CHKOP, Alloc>::mod_type::first_type &
        m ( y.mod().first );

        if ( m != Rational<mpz_class, GCD, CHKOP, Alloc>::zero_ &&
                mpz_perfect_square_p ( m.get_mpz_t() ) ) {

            typename Rational<mpz_class, GCD, CHKOP, Alloc>::integer_type psq;
            mpz_sqrt ( psq.get_mpz_t(), m.get_mpz_t() );
            return Rational<mpz_class, GCD, CHKOP, Alloc> ( psq, Rational<mpz_class, GCD, CHKOP,
                    Alloc>::one_ );
        }

        return x;
    }
};

/**
 * @ingroup gmp
 * @brief Rational class based on the GMP library
 */
typedef Rational<mpz_class, GCD_gmp, NO_OPERATOR_CHECK> gmp_rational;

#ifndef GMP_HERON_DIGITS
/**
 * @ingroup gmp
 * @def GMP_HERON_DIGITS
 *
 * @brief Upper bound in digits of the denominator of square root approximations
 *
 * @see Commons::Math::Rational::sqrt()
 */
#define GMP_HERON_DIGITS 28u
#endif

template<> inline bool SQRT_HERON_ITERATE<gmp_rational>::operator() ( const gmp_rational &p,
        const mpz_class &, const mpz_class & ) const {

    const gmp_rational::mod_type &m ( p.mod() );

    return ! ( gmp_rational::isInteger ( m ) && mpz_perfect_square_p ( m.first.get_mpz_t() ) );
}

template<> inline bool SQRT_HERON_ITERATE<gmp_rational>::operator() ( const gmp_rational &x,
        const gmp_rational & ) const {
    return mpz_sizeinbase ( x.denominator().get_mpz_t(), 10 ) < GMP_HERON_DIGITS;
}

template<> struct CFRationalTraits<mpz_class> {
    typedef gmp_rational rational_type;
};

template<template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP, template<typename> class Alloc>
struct _remquo<mpz_class, GCD, CHKOP, Alloc> {

    mpz_class operator() ( const mpz_class &x, const mpz_class &y, mpz_class &quo ) const {
        mpz_class r;
        mpz_fdiv_qr ( quo.get_mpz_t(), r.get_mpz_t(), x.get_mpz_t(), y.get_mpz_t() );
        return r;
    }
};

}

}

#ifdef HAVE_MPREAL_H
inline mpz_class ceil ( const mpf_class &f ) {
    mpf_class aux;
    mpf_ceil ( aux.get_mpf_t(), f.get_mpf_t() );
    return mpz_class ( aux );
}

inline mpf_class log10 ( const mpz_class &z ) {
    mpf_class f;
    mpfr_get_f ( f.get_mpf_t(), mpfr::log10 ( mpfr::mpreal ( z.get_mpz_t() ) ).mpfr_ptr(),
                 mpfr::mpreal::get_default_rnd() );
    return f;
}

inline mpz_class pow10 ( const mpf_class &f ) {
    mpz_class z;
    mpfr_get_z ( z.get_mpz_t(), mpfr::pow ( 10.0, mpfr::mpreal ( f.get_mpf_t() ) ).mpfr_ptr(),
                 mpfr::mpreal::get_default_rnd() );
    return z;
}

inline mpz_class floor ( const mpz_class &z ) {
    mpz_class r;
    mpfr::mpreal rf, f ( z.get_mpz_t() );
    mpfr_floor ( rf.mpfr_ptr(), f.mpfr_ptr() );
    mpfr_get_z ( r.get_mpz_t(), rf.mpfr_ptr(), mpfr::mpreal::get_default_rnd() );
    return r;
}
#else
inline mpz_class ceil ( const mpf_class &f ) {
    return std::ceil ( f.get_d() );
}

inline mpf_class log10 ( const mpz_class &z ) {
    return ( std::log10 ( z.get_d() ) );
}

inline mpz_class pow10 ( const mpf_class &f ) {
    return std::pow ( 10, f.get_d() );
}

inline mpz_class floor ( const mpz_class &z ) {
    mpf_class r, f ( z );
    mpf_floor ( r.get_mpf_t(), f.get_mpf_t() );
    mpz_class aux;
    mpz_set_f ( aux.get_mpz_t(), r.get_mpf_t() );
    return aux;
}
#endif

#endif /* COMMONS_MATH_GMP_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
