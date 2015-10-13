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
 * @author Heiko Schäfer <heiko@rangun.de>
 * @copyright 2015 by Heiko Schäfer <heiko@rangun.de>
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
 *
 * @warning You cannot use @ref cln and @ref gmp in the same compilation unit
 */

#ifndef COMMONS_MATH_CLN_RATIONAL_H
#define COMMONS_MATH_CLN_RATIONAL_H

#ifdef COMMONS_MATH_GMP_RATIONAL_H
#error "cln_rational.h and gmp_rational.h cannot coexist in the same compilation unit"
#endif

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
 * The default precision suffix
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
 * The @c EPSILON used for approximating a float
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

    static cln::cl_I max() {
        return cln::cl_I();
    }

    static cln::cl_I min() {
        return cln::cl_I();
    }

    static cln::cl_I epsilon() {
        return cln::cl_I();
    }
};

#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
template<> struct divides<cln::cl_I> :
    public binary_function<cln::cl_I, cln::cl_I, cln::cl_I> {
    inline result_type operator() ( const first_argument_type &x,
                                    const second_argument_type &y ) const {
        return cln::truncate1 ( x, y );
    }
};

template<> struct modulus<cln::cl_I> :
    public binary_function<cln::cl_I, cln::cl_I, cln::cl_I> {
    inline result_type operator() ( const first_argument_type &x,
                                    const second_argument_type &y ) const {
        return cln::truncate2 ( x, y ).remainder;
    }
};
#pragma GCC diagnostic pop

}

namespace Commons {

namespace Math {

template<> struct TYPE_CONVERT<float> {

    inline explicit TYPE_CONVERT ( float v ) : val ( v ) {}

    template<class U>
    inline U convert() const {
        return cln::floor1 ( val );
    }

private:
    float val;
};

template<> struct TYPE_CONVERT<double> {

    inline explicit TYPE_CONVERT ( double v ) : val ( v ) {}

    template<class U>
    inline U convert() const {
        return cln::floor1 ( val );
    }

private:
    double val;
};

template<> struct TYPE_CONVERT<cln::cl_F> {

    inline explicit TYPE_CONVERT ( const cln::cl_F &v ) : val ( v ) {}

    inline explicit TYPE_CONVERT ( const cln::cl_I &v ) : val ( cln::double_approx ( v ) ) {}

    template<class U>
    inline U convert() const {
        return cln::floor1 ( val );
    }

private:
    const cln::cl_F val;
};

template<> struct TYPE_CONVERT<cln::cl_I> {

    inline explicit TYPE_CONVERT ( const cln::cl_I &v ) : val ( v ) {}

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

template<> inline cln::cl_I TYPE_CONVERT<cln::cl_I>::convert<cln::cl_I>() const {
    return val;
}

template<> struct EPSILON<cln::cl_F> {
    inline static const cln::cl_F value() {
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

    inline cln::cl_I operator() ( const cln::cl_I &a, const cln::cl_I &b ) const {
        return GCD_cln<cln::cl_I, true, CHKOP, CONV>() ( a, b );
    }

};

template<template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_cln<cln::cl_I, true, CHKOP, CONV> {

    inline cln::cl_I operator() ( const cln::cl_I &a, const cln::cl_I &b ) const {
        return cln::gcd ( a, b );
    }

};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _lcm<cln::cl_I, GCD, CHKOP, true> {
    inline cln::cl_I operator() ( const cln::cl_I &a, const cln::cl_I &b ) const {
        return cln::lcm ( a, b );
    }
};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _lcm<cln::cl_I, GCD, CHKOP, false> {

    inline cln::cl_I operator() ( const cln::cl_I &a, const cln::cl_I &b ) const {
        return _lcm<cln::cl_I, GCD, CHKOP, true>() ( a, b );
    }
};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _abs<cln::cl_I, GCD, CHKOP, true> {

    inline Rational<cln::cl_I, GCD, CHKOP> operator()
    ( const Rational<cln::cl_I, GCD, CHKOP> &r ) const {
        return Rational<cln::cl_I, GCD, CHKOP> ( cln::abs ( r.numerator() ), r.denominator() );
    }
};

template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _abs<cln::cl_I, GCD, CHKOP, false> {

    inline Rational<cln::cl_I, GCD, CHKOP> operator()
    ( const Rational<cln::cl_I, GCD, CHKOP> &r ) const {
        return _abs<cln::cl_I, GCD, CHKOP, true>() ( r );
    }
};

template<>
inline cln::cl_F _approxFract<cln::cl_I, GCD_cln, NO_OPERATOR_CHECK, cln::cl_F, true,
EPSILON, TYPE_CONVERT>::abs ( const cln::cl_F &nt ) const {
    return cln::abs ( nt );
}

/**
 * @ingroup cln
 * @brief Rational class based on the CLN library
 */

typedef Rational<cln::cl_I, Commons::Math::GCD_cln, Commons::Math::NO_OPERATOR_CHECK> cln_rational;

}

}

namespace cln {

inline cl_I floor ( const cl_F &f ) {
    return floor1 ( f );
}

}

#endif /* COMMONS_MATH_CLN_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
