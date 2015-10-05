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

#ifndef GMPTESTCASE_H
#define GMPTESTCASE_H

#include <gmpxx.h>

#include <cppunit/extensions/HelperMacros.h>

#include "rational.h"

namespace std {

#if (__GNU_MP_VERSION * 10000 + __GNU_MP_VERSION_MINOR * 100 + __GNU_MP_VERSION_PATCHLEVEL) \
          < 50100

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

template<>
struct numeric_limits<mpf_class> {
    static const bool is_specialized = true;
    static mpf_class min() {
        return mpf_class();
    }
    static mpf_class max() {
        return mpf_class();
    }
    static mpf_class lowest() {
        return mpf_class();
    }
    static const int digits = 0;
    static const int digits10 = 0;
    static const int max_digits10 = 0;
    static const bool is_signed = true;
    static const bool is_integer = false;
    static const bool is_exact = false;
    static const int radix = 2;
    static mpf_class epsilon() {
        return mpf_class ();
    }
    static mpf_class round_error() {
        return mpf_class();
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
    static mpf_class infinity() {
        return mpf_class();
    }
    static mpf_class quiet_NaN() {
        return mpf_class();
    }
    static mpf_class signaling_NaN() {
        return mpf_class();
    }
    static mpf_class denorm_min() {
        return mpf_class();
    }
    static const bool is_iec559 = false;
    static const bool is_bounded = false;
    static const bool is_modulo = false;
    static const bool traps = false;
    static const bool tinyness_before = false;
    static const float_round_style round_style = round_indeterminate;
};

#endif
}

namespace Commons {

namespace Math {

template<> struct TYPE_CONVERT<mpz_class> {

    inline explicit TYPE_CONVERT ( const mpz_class &v ) : val ( v ) {}

    inline operator mpf_class() const {
        return mpf_class ( val );
    }

private:
    const mpz_class &val;
};

template<> mpf_class EPSILON<mpf_class>::value() {
    return mpf_class ( "1e-21", 30, 10 );
}

template<>
inline mpf_class _approxFract<mpz_class, GCD_euclid, NO_OPERATOR_CHECK, mpf_class, true,
EPSILON, TYPE_CONVERT>::abs ( const mpf_class &nt ) const {
    return ::abs ( nt );
}

}

}

#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
class GMPTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( GMPTest );
    CPPUNIT_TEST ( testConstruct );
    CPPUNIT_TEST ( testConstructFrom_mpf_class );
    CPPUNIT_TEST ( testAlgorithm );
    CPPUNIT_TEST_SUITE_END();

public:
    typedef mpz_class rational_type;

    GMPTest();

    void setUp();
    void tearDown();

    void testConstruct();
    void testConstructFrom_mpf_class();
    void testAlgorithm();

private:
    typedef Commons::Math::Rational<rational_type, Commons::Math::GCD_euclid> unchecked_sqrt;
    typedef std::vector<unchecked_sqrt> rat_vector_sqrt;
    typedef std::vector<Commons::Math::Rational<rational_type> > rat_vector;

    unchecked_sqrt m_sqrt2;
    rat_vector_sqrt m_twosqrt;
    rat_vector m_onethird;
    rat_vector m_oneseventh;
};
#pragma GCC diagnostic pop

#endif /* GMPTESTCASE_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
