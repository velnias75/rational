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
 * @author Heiko Schäfer <heiko@rangun.de>
 * @copyright 2015-2016 by Heiko Schäfer <heiko@rangun.de>
 *
 * This header contains specialization for use of the
 * [CLN - Class Library for Numbers](http://www.ginac.de/CLN/) as underlying
 * storage type.
 *
 * @defgroup cln CLN - Class Library for Numbers extensions
 *
 * The header `cln_rational.h` contains specializations especially for the
 * [CLN - Class Library for Numbers](http://www.ginac.de/CLN/)
 * as underlying storage type.\n
 * \n If you use the *CLN extensions*, you'll need to link your
 * application with `-lcln`
 */

#ifndef COMMONS_MATH_CLN_RATIONAL_H
#define COMMONS_MATH_CLN_RATIONAL_H

#include <sstream>

#include <cln/integer_io.h>
#include <cln/modinteger.h>
#include <cln/integer.h>
#include <cln/float.h>

#include "rational.h"

#ifndef CLN_PRECISION
/**
 * @ingroup cln
 * @def CLN_PRECISION
 *
 * @brief The default precision suffix
 *
 * @see CLN_EPSILON
 * @see Commons::Math::EPSILON
 * @see Commons::Math::TYPE_CONVERT
 */
#define CLN_PRECISION "30"
#endif


#ifndef CLN_EPSILON
/**
 * @ingroup cln
 * @def CLN_EPSILON
 *
 * @brief The @c EPSILON used for approximating a float
 *
 * @see CLN_PRECISION
 * @see Commons::Math::EPSILON
 *
 * This define is passed @em as @em is to the contructor of @c cln::cl_F
 */
#define CLN_EPSILON "1L-16_" CLN_PRECISION
#endif

namespace std {

template<> struct numeric_limits<cln::cl_I> {

    static bool const is_specialized = true;
    static bool const is_signed = true;
    static bool const is_integer = true;
    static bool const is_exact = true;

    static const cln::cl_I max() {
        return cln::cl_I();
    }

    static const cln::cl_I min() {
        return cln::cl_I();
    }

    static const cln::cl_I epsilon() {
        return cln::cl_I();
    }
};

#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
template<> struct divides<cln::cl_I> :
    public binary_function<cln::cl_I, cln::cl_I, cln::cl_I> {
    result_type operator() ( const first_argument_type &x,
                             const second_argument_type &y ) const {
        return cln::truncate1 ( x, y );
    }
};

template<> struct modulus<cln::cl_I> :
    public binary_function<cln::cl_I, cln::cl_I, cln::cl_I> {
    result_type operator() ( const first_argument_type &x,
                             const second_argument_type &y ) const {
        return cln::truncate2 ( x, y ).remainder;
    }
};
#pragma GCC diagnostic pop

}

namespace Commons {

namespace Math {

template<> struct ExpressionEvalTraits<cln::cl_I> {
    typedef cln::cl_F NumberType;
};

template<> struct DecomposeBaseTraits<cln::cl_I, true> {
    typedef cln::cl_I digit_type;
    enum { Base = 10 };
};

template<> struct _type_round_helper<cln::cl_I> {
    cln::cl_I operator() ( const cln::cl_I &tr ) const {
        return tr;
    }
};

template<> inline cln::cl_I TYPE_CONVERT<long double>::convert<cln::cl_I>() const {
    std::ostringstream os;
    os.precision ( std::numeric_limits<long double>::digits );
    os << val;
    return os.str().c_str();
}

template<> inline cln::cl_F TYPE_CONVERT<std::string>::convert<cln::cl_F>() const {
    return std::string ( val ).append ( "L0_" ).append ( CLN_PRECISION ).c_str();
}

template<> inline cln::cl_F TYPE_CONVERT<const char *>::convert<cln::cl_F>() const {
    return ( len ? std::string ( val, len ) : std::string ( val ) ).append ( "L0_" ).
           append ( CLN_PRECISION ).c_str();
}

template<> struct TYPE_CONVERT<cln::cl_F> {

    RATIONAL_NOCOPYASSIGN ( TYPE_CONVERT<cln::cl_F> );

    explicit TYPE_CONVERT ( const cln::cl_F &v ) : val ( v ) {}

    explicit TYPE_CONVERT ( const cln::cl_I &v ) : val ( cln::double_approx ( v ) ) {}

    template<class U> U convert() const {
        return cln::floor1 ( val );
    }

private:
    const cln::cl_F val;
};

template<> inline cln::cl_F TYPE_CONVERT<cln::cl_F>::convert<cln::cl_F>() const {
    return val;
}

template<> struct TYPE_CONVERT<cln::cl_I> {

    RATIONAL_NOCOPYASSIGN ( TYPE_CONVERT<cln::cl_I> );

    explicit TYPE_CONVERT ( const cln::cl_I &v ) : val ( v ) {}

    template<class U> U convert() const;

private:
    const cln::cl_I &val;
};

template<class U>
U TYPE_CONVERT<cln::cl_I>::convert() const {
    std::ostringstream os;
    os << val << "L+0_" << CLN_PRECISION;
    return os.str().c_str();
}

template<> inline float TYPE_CONVERT<cln::cl_I>::convert<float>() const {
    return cln::float_approx ( val );
}

template<> inline double TYPE_CONVERT<cln::cl_I>::convert<double>() const {
    return cln::double_approx ( val );
}

template<> inline long double TYPE_CONVERT<cln::cl_I>::convert<long double>() const {
    return static_cast<long double> ( cln::double_approx ( val ) );
}

template<> inline cln::cl_I TYPE_CONVERT<cln::cl_I>::convert<cln::cl_I>() const {
    return val;
}

template<> inline cln::cl_I TYPE_CONVERT<float>::convert<cln::cl_I>() const {
    return cln::floor1 ( val );
}

template<> inline cln::cl_I TYPE_CONVERT<double>::convert<cln::cl_I>() const {
    return cln::floor1 ( val );
}

template<> struct EPSILON<cln::cl_F> {

    RATIONAL_NOCOPYASSIGN ( EPSILON<cln::cl_F> );

    static const cln::cl_F value() {
        return CLN_EPSILON;
    }
};

/**
 * @ingroup cln
 * @ingroup gcd
 * @brief CLN GCD algorithm
 *
 * The gcd algorithm implemented in the CLN library
 *
 * @note this is currently the @b only supported GCD for use with @ref cln
 *
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, bool IsSigned, template<class, typename = T, bool = IsSigned> class CHKOP,
         template<typename = T> class CONV> struct GCD_cln;

template<template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_cln<cln::cl_I, false, CHKOP, CONV> {

    cln::cl_I operator() ( const cln::cl_I &a, const cln::cl_I &b ) const {
        return GCD_cln<cln::cl_I, true, CHKOP, CONV>() ( a, b );
    }

};

template<template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_cln<cln::cl_I, true, CHKOP, CONV> {

    cln::cl_I operator() ( const cln::cl_I &a, const cln::cl_I &b ) const {
        return cln::gcd ( a, b );
    }

};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc> struct _lcm<cln::cl_I, GCD, CHKOP, Alloc, true> {

    cln::cl_I operator() ( const cln::cl_I &a, const cln::cl_I &b ) const {
        return cln::lcm ( a, b );
    }
};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc> struct _lcm<cln::cl_I, GCD, CHKOP, Alloc, false> {

    cln::cl_I operator() ( const cln::cl_I &a, const cln::cl_I &b ) const {
        return _lcm<cln::cl_I, GCD, CHKOP, Alloc, true>() ( a, b );
    }
};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc> struct _abs<cln::cl_I, GCD, CHKOP, Alloc, true> {

    Rational<cln::cl_I, GCD, CHKOP, Alloc> operator()
    ( const Rational<cln::cl_I, GCD, CHKOP, Alloc> &r ) const {
        return Rational<cln::cl_I, GCD, CHKOP, Alloc> ( cln::abs ( r.numerator() ),
                r.denominator() );
    }
};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc> struct _abs<cln::cl_I, GCD, CHKOP, Alloc, false> {

    Rational<cln::cl_I, GCD, CHKOP, Alloc> operator()
    ( const Rational<cln::cl_I, GCD, CHKOP, Alloc> &r ) const {
        return _abs<cln::cl_I, GCD, CHKOP, Alloc, true>() ( r );
    }
};

template<template<typename> class EPSILON> struct _approxUtils<cln::cl_F, EPSILON> {

    static bool approximated ( const cln::cl_F &af, const cln::cl_F &nt ) {
        return abs ( af - nt ) < eps_;
    }

    static cln::cl_F reciprocal ( const cln::cl_F &x ) {
        return one_ / x;
    }

    const static cln::cl_F eps_;

private:
    static cln::cl_F abs ( const cln::cl_F &nt ) {
        return cln::abs ( nt );
    }

private:
    const static cln::cl_F one_;
};

template<template<typename> class EPSILON>
const cln::cl_F _approxUtils<cln::cl_F, EPSILON>::eps_ ( EPSILON<cln::cl_F>::value() );

template<template<typename> class EPSILON>
const cln::cl_F _approxUtils<cln::cl_F, EPSILON>::one_ ( 1.0 );

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc> struct _swapSign<cln::cl_I, GCD, CHKOP, Alloc, true> {

    Rational<cln::cl_I, GCD, CHKOP, Alloc> &
    operator() ( Rational<cln::cl_I, GCD, CHKOP, Alloc> &r ) const {

        if ( cln::signum ( r.m_denom ) < zero_ ) {
            r.m_numer = -r.m_numer;
            r.m_denom = -r.m_denom;
        }

        return r;
    }

private:
    static const cln::cl_I zero_;
};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc>
const cln::cl_I _swapSign<cln::cl_I, GCD, CHKOP, Alloc, true>::zero_ ( 0 );

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         template<typename> class Alloc> struct _psq<cln::cl_I, GCD, CHKOP, Alloc> {

    Rational<cln::cl_I, GCD, CHKOP, Alloc> operator()
    ( const Rational<cln::cl_I, GCD, CHKOP, Alloc> &x,
      const Rational<cln::cl_I, GCD, CHKOP, Alloc> &y ) const {

        const typename Rational<cln::cl_I, GCD, CHKOP, Alloc>::mod_type::first_type &
        m ( y.mod().first );

        typename Rational<cln::cl_I, GCD, CHKOP, Alloc>::integer_type psq;

        if ( m != Rational<cln::cl_I, GCD, CHKOP, Alloc>::zero_ && cln::sqrtp ( m, &psq ) ) {
            return Rational<cln::cl_I, GCD, CHKOP, Alloc> ( psq, Rational<cln::cl_I, GCD, CHKOP,
                    Alloc>::one_ );
        }

        return x;
    }
};

/**
 * @ingroup cln
 * @brief Rational class based on the CLN library
 */

typedef Rational<cln::cl_I, Commons::Math::GCD_cln, Commons::Math::NO_OPERATOR_CHECK> cln_rational;

#ifndef CLN_HERON_DIGITS
/**
 * @ingroup cln
 * @def CLN_HERON_DIGITS
 *
 * @brief Upper bound in digits of the denominator of square root approximations
 *
 * @see Commons::Math::Rational::sqrt()
 */
#define CLN_HERON_DIGITS 28u
#endif

template<> inline bool SQRT_HERON_ITERATE<cln_rational>::operator() ( const cln_rational &p,
        const cln::cl_I &, const cln::cl_I & ) const {

    const cln_rational::mod_type &m ( p.mod() );
    typename cln_rational::integer_type psq;

    return ! ( cln_rational::isInteger ( m ) && cln::sqrtp ( m.first, &psq ) );
}

template<> inline bool SQRT_HERON_ITERATE<cln_rational>::operator() ( const cln_rational &x,
        const cln_rational & ) const {
    return cln::integer_length ( x.denominator() ) < ( CLN_HERON_DIGITS * 3u );
}

template<> struct CFRationalTraits<cln::cl_I> {
    typedef cln_rational rational_type;
};

template<template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP, template<typename> class Alloc>
struct _remquo<cln::cl_I, GCD, CHKOP, Alloc> {

    cln::cl_I operator() ( const cln::cl_I &x, const cln::cl_I &y, cln::cl_I &quo ) const {

        const cln::cl_I_div_t &d ( cln::floor2 ( x, y ) );

        quo = d.quotient;
        return d.remainder;
    }
};

}

}

namespace cln {

inline cl_I floor ( const cl_I &i ) {
    return floor1 ( cln::double_approx ( i ) );
}

inline cl_I floor ( const cl_F &f ) {
    return floor1 ( f );
}

inline cl_I ceil ( const cl_F &f ) {
    return ceiling1 ( f );
}

inline cl_F log10 ( const cl_I &i ) {
    return std::log10 ( cln::double_approx ( i ) );
}

inline cl_I pow10 ( const cl_I &i ) {
    return static_cast<unsigned long> ( std::pow ( 10, cln::double_approx ( i ) ) );
}

}

#endif /* COMMONS_MATH_CLN_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
