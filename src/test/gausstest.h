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

#ifndef GAUSSTESTCASE_H
#define GAUSSTESTCASE_H

#include <complex>

#include <cppunit/extensions/HelperMacros.h>

#include "rational.h"

namespace std {

template<> struct numeric_limits<complex<long> > {

    static bool const is_specialized = true;
    static bool const is_signed = true;
    static bool const is_integer = true;
    static bool const is_exact = true;

    static const complex<long> max() {
        return complex<long>();
    }

    static const complex<long> min() {
        return complex<long>();
    }

    static const complex<long> epsilon() {
        return complex<long>();
    }
};

#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
template<> struct modulus<complex<long> > : binary_function<complex<long>, complex<long>,
        complex<long> > {

    inline result_type operator() ( const first_argument_type &x,
                                    const second_argument_type &y ) const {

        long q1 = x.real() / y.real();
        long r1 = x.real() - ( q1 * y.real() );

        long q2 = x.imag() / y.imag();
        long r2 = x.imag() - ( q2 * y.imag() );

        return result_type ( r1, r2 );
    }
};
#pragma GCC diagnostic pop

}

namespace Commons {

namespace Math {

template<template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_euclid<std::complex<long>, true, CHKOP, CONV> {

    std::complex<long> operator() ( const std::complex<long> &a,
                                    const std::complex<long> &b ) const;

private:
    static const std::complex<long> zero_;
};

template<template<class, typename, bool> class CHKOP, template<typename> class CONV>
const std::complex<long> GCD_euclid<std::complex<long>, true, CHKOP, CONV>::zero_
    = std::complex<long>();

template<template<class, typename, bool> class CHKOP, template<typename> class CONV>
std::complex<long> GCD_euclid<std::complex<long>, true, CHKOP, CONV>::operator()
( const std::complex<long> &a, const std::complex<long> &b ) const {
    return std::abs ( GCD_euclid<std::complex<long>, false, CHKOP, CONV>() ( a, b ) );
}

#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic push
template<template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _swapSign<std::complex<long>, GCD, CHKOP, true> {

    inline Rational<std::complex<long>, GCD, CHKOP> &operator() ( Rational<std::complex<long>, GCD,
            CHKOP> &r ) const {

        if ( std::abs ( r.denominator() != r.denominator() ) ) {
            r.m_numer = -r.m_numer;
            r.m_denom = -r.m_denom;
        }

        return r;
    }
};
#pragma GCC diagnostic pop

}

}

#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
class GaussTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( GaussTest );
    CPPUNIT_TEST ( testConstruct );
    CPPUNIT_TEST_SUITE_END();

public:
    typedef Commons::Math::Rational<std::complex<long>, Commons::Math::GCD_euclid> gauss_rational;

    GaussTest();

    void setUp();
    void tearDown();

    void testConstruct();
};

#endif /* GAUSSTESTCASE_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
