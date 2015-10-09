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

#include <numeric>

#include "gmptest.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( GMPTest );

using namespace Commons::Math;

GMPTest::GMPTest() : CppUnit::TestFixture(), m_sqrt2(), m_twosqrt(), m_onethird(),
    m_oneseventh() {}

void GMPTest::setUp() {

    m_sqrt2 = unchecked_sqrt ( mpf_class ( std::sqrt ( mpf_class ( 2.0 ) ) ) );

    std::fill_n ( std::back_inserter ( m_twosqrt ), 2, m_sqrt2 );
    std::fill_n ( std::back_inserter ( m_onethird ), 3, rat_vector::value_type ( 1, 3 ) );
    std::fill_n ( std::back_inserter ( m_oneseventh ), 7, rat_vector::value_type ( 1, 7 ) );
}

void GMPTest::tearDown() {}

void GMPTest::testConstruct() {

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( gmp_rational r ( 1, 0 ), std::domain_error );
#endif

    CPPUNIT_ASSERT_EQUAL (
        std::string ( "4016738247910921445183770562792527230871790482730798620920" \
                      "3242394032674573345714107911694514536898289177933847227109" \
                      "12324725149799230371316644957" ), m_sqrt2.numerator().get_str() );

    CPPUNIT_ASSERT_EQUAL (
        std::string ( "2840262853349184214224624055939823380440794351393149143229" \
                      "9345538280058115643575803799519938893676390738702021028792" \
                      "52515818756635810738450587537" ), m_sqrt2.denominator().get_str() );

    const mpf_class &a ( gmp_rational ( 1, 2 ) );
    const mpf_class &b ( gmp_rational ( 1, -2 ) );
    const mpf_class &c ( gmp_rational ( -1, 2 ) );
    const mpf_class &d ( gmp_rational ( -1, -2 ) );

    CPPUNIT_ASSERT_EQUAL ( 0.5, a.get_d() );
    CPPUNIT_ASSERT_EQUAL ( -0.5, b.get_d() );
    CPPUNIT_ASSERT_EQUAL ( -0.5, c.get_d() );
    CPPUNIT_ASSERT_EQUAL ( 0.5, d.get_d() );

    CPPUNIT_ASSERT_EQUAL ( -3l, gmp_rational ( 6, -8 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4l, gmp_rational ( 6, -8 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 7l, gmp_rational ( 14, 24 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, gmp_rational ( 14, 24 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 7l, gmp_rational ( 2, 1, 3 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 3l, gmp_rational ( 2, 1, 3 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 86l, gmp_rational ( 18, 4, -5 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 5l, gmp_rational ( 18, 4, -5 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 86l, gmp_rational ( 18, -4, 5 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 5l, gmp_rational ( 18, -4, 5 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( -86l, gmp_rational ( -18, 4, 5 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 5l, gmp_rational ( -18, 4, 5 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( -94l, gmp_rational ( -18, 4, -5 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 5l, gmp_rational ( -18, 4, -5 ).denominator().get_si() );
}

void GMPTest::testConstructFromDouble() {

    const Rational<rational_type, GCD_euclid> &p ( 19.0 / 51.0 );

    CPPUNIT_ASSERT_EQUAL ( 19l, p.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 51l, p.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &q ( 516901.0 / 740785.0 );

    CPPUNIT_ASSERT_EQUAL ( 516901l, q.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 740785l, q.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &r ( -0.7391304347826086 );

    CPPUNIT_ASSERT_EQUAL ( -17l, r.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 23l, r.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &s ( 0.0 );

    CPPUNIT_ASSERT_EQUAL ( 0l, s.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, s.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &pi ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( 245850922l, pi.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 78256779l, pi.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &t ( 1.0 );

    CPPUNIT_ASSERT_EQUAL ( 1l, t.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, t.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &u ( 2.0 );

    CPPUNIT_ASSERT_EQUAL ( 2l, u.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, u.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &v ( -8 );

    CPPUNIT_ASSERT_EQUAL ( -8l, v.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, v.denominator().get_si() );
}

void GMPTest::testConstructFrom_mpf_class() {

    const Rational<rational_type, GCD_euclid> &p ( mpf_class ( mpf_class ( 19.0 ) /
            mpf_class ( 51.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 19l, p.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 51l, p.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &q ( mpf_class ( mpf_class ( 516901.0 ) /
            mpf_class ( 740785.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 516901l, q.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 740785l, q.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &r ( mpf_class ( -0.7391304347826086, 65 ) );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-665749510133023" ), r.numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "900719925474090" ), r.denominator().get_str() );

    const Rational<rational_type, GCD_euclid> &s ( mpf_class ( 0.0 ) );

    CPPUNIT_ASSERT_EQUAL ( 0l, s.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, s.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &pi ( mpf_class ( M_PI ) );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "9978066541" ), pi.numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "3176117225" ), pi.denominator().get_str() );

    const Rational<rational_type, GCD_euclid> &t ( mpf_class ( 1.0 ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, t.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, t.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &u ( mpf_class ( 2.0 ) );

    CPPUNIT_ASSERT_EQUAL ( 2l, u.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, u.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> &v ( mpf_class ( -8 ) );

    CPPUNIT_ASSERT_EQUAL ( -8l, v.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, v.denominator().get_si() );
}

void GMPTest::testAddition() {

    const Rational<rational_type, GCD_euclid_fast> a ( 17, 21 );
    const Rational<rational_type, GCD_euclid_fast> b ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( 31l, ( a + b ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 15l, ( a + b ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 31l, ( b + a ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 15l, ( b + a ).denominator().get_si() );

    const Rational<rational_type, GCD_euclid> c ( 1, 6 );
    const Rational<rational_type, GCD_euclid> d ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( 3l, ( c + d ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 10l, ( c + d ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 3l, ( d + c ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 10l, ( d + c ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 2l, ( +d ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 15l, ( +d ).denominator().get_si() );

    const Rational<rational_type, GCD_stein> knuth_a ( 7, 66 );
    const Rational<rational_type, GCD_stein> knuth_b ( 17, 12 );

    CPPUNIT_ASSERT_EQUAL ( 67l, ( knuth_a + knuth_b ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 44l, ( knuth_a + knuth_b ).denominator().get_si() );
}

void GMPTest::testSubtraction() {

    const Rational<rational_type, GCD_euclid> a ( 17, 21 );
    const Rational<rational_type, GCD_euclid> b ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( -47l, ( a - b ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 105l, ( a - b ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 0l, ( a - a ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, ( a - a ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 47l, ( b - a ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 105l, ( b - a ).denominator().get_si() );

    const Rational<rational_type, GCD_euclid> c ( 1, 6 );
    const Rational<rational_type, GCD_euclid> d ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( 1l, ( c - d ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 30l, ( c - d ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( -1l, ( d - c ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 30l, ( d - c ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( -2l, ( -d ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 15l, ( -d ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 2l, ( d ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 15l, ( d ).denominator().get_si() );
}

void GMPTest::testMultiplication() {

    const gmp_rational a ( 2, 8 );
    const gmp_rational b ( 7, 3 );

    CPPUNIT_ASSERT_EQUAL ( 7l, ( a * b ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, ( a * b ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 7l, ( b * a ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, ( b * a ).denominator().get_si() );
}

void GMPTest::testInvert() {

    CPPUNIT_ASSERT_EQUAL ( 7l, gmp_rational ( 161, 49 ).invert().numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 23l,
                           gmp_rational ( 161, 49 ).invert().denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 7l, gmp_rational ( 161, 49 ).inverse().numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 23l,
                           gmp_rational ( 161, 49 ).inverse().denominator().get_si() );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( gmp_rational ().invert(), std::domain_error );
    CPPUNIT_ASSERT_THROW ( gmp_rational ().inverse(), std::domain_error );
#endif
}

void GMPTest::testDivision() {

    const gmp_rational a ( 2, 8 );
    const gmp_rational b ( 7, 3 );
    const gmp_rational c ( 0, 1 );
    const gmp_rational d ( -7, -3 );

    CPPUNIT_ASSERT_EQUAL ( 3l, ( a / b ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 28l, ( a / b ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 28l, ( b / a ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 3l, ( b / a ).denominator().get_si() );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( a / c, std::domain_error );
    CPPUNIT_ASSERT_THROW ( a / ( b - d ), std::domain_error );
#endif
}

void GMPTest::testModulo() {

    gmp_rational a ( 8, 1 );

    a %= gmp_rational ( 3, 1 );

    CPPUNIT_ASSERT_EQUAL ( 2l, ( a ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, ( a ).denominator().get_si() );

    gmp_rational c ( 41, 7 );
    c %= gmp_rational ( 3, 2 );

    CPPUNIT_ASSERT_EQUAL ( 19l, ( c ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 14l, ( c ).denominator().get_si() );

    gmp_rational d ( 542, 84 );
    gmp_rational e ( -65, 28 );

    CPPUNIT_ASSERT_EQUAL ( -43l, ( d % e ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 84l, ( d % e ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 347l, ( e % d ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 84l, ( e % d ).denominator().get_si() );

    gmp_rational h ( 11, 4 );

    CPPUNIT_ASSERT_EQUAL ( 2l, h.mod().first.get_si() );
    CPPUNIT_ASSERT_EQUAL ( 3l, h.mod().second.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4l, h.mod().second.denominator().get_si() );

    gmp_rational i ( 11, -4 );

    CPPUNIT_ASSERT_EQUAL ( -2l, i.mod().first.get_si() );
    CPPUNIT_ASSERT_EQUAL ( 3l, i.mod().second.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4l, i.mod().second.denominator().get_si() );

    gmp_rational j ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( 2l, j.mod().first.get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, j.mod().second.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4l, j.mod().second.denominator().get_si() );

    gmp_rational k ( -18, 8 );

    CPPUNIT_ASSERT_EQUAL ( -2l, k.mod().first.get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, k.mod().second.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4l, k.mod().second.denominator().get_si() );

    gmp_rational l ( 1, 8 );

    CPPUNIT_ASSERT_EQUAL ( 0l, l.mod().first.get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, l.mod().second.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 8l, l.mod().second.denominator().get_si() );
}

void GMPTest::testIncDec() {

    gmp_rational a ( 2, 4 );

    CPPUNIT_ASSERT_EQUAL ( 3l, ( ++a ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 2l, ( a++ ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 5l, a.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 2l, a.denominator().get_si() );

    gmp_rational b ( 2, 4 );

    CPPUNIT_ASSERT_EQUAL ( -1l, ( --b ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 2l, ( b-- ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( -3l, b.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 2l, b.denominator().get_si() );
}

void GMPTest::testRelOps() {

    const gmp_rational a ( 1, 4 );
    const gmp_rational b ( 1, 2 );

    CPPUNIT_ASSERT ( a < b );
    CPPUNIT_ASSERT ( a <= b );

    CPPUNIT_ASSERT ( b > a );
    CPPUNIT_ASSERT ( b >= a );

    const gmp_rational c ( 2, 4 );

    CPPUNIT_ASSERT ( c == b );
    CPPUNIT_ASSERT ( b == c );

    CPPUNIT_ASSERT ( a != b );
    CPPUNIT_ASSERT ( b != a );

    CPPUNIT_ASSERT ( b <= c );
    CPPUNIT_ASSERT ( c <= b );
    CPPUNIT_ASSERT ( b >= c );
    CPPUNIT_ASSERT ( c >= b );

    const gmp_rational d ( 2, 4 );
    const gmp_rational e ( 2, -4 );

    CPPUNIT_ASSERT ( d > e );
    CPPUNIT_ASSERT ( e < d );

    const gmp_rational f ( -2, 4 );

    CPPUNIT_ASSERT ( f == e );
    CPPUNIT_ASSERT ( f >= e );
    CPPUNIT_ASSERT ( f <= e );

    CPPUNIT_ASSERT ( e == f );
    CPPUNIT_ASSERT ( e >= f );
    CPPUNIT_ASSERT ( e <= f );

    const gmp_rational g ( -3, 4 );

    CPPUNIT_ASSERT ( g < d );
    CPPUNIT_ASSERT ( d > g );
}

void GMPTest::testString() {

    const gmp_rational h ( 11, 4 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/4" ), h.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "2 3/4" ), h.str ( true ) );

    const gmp_rational i ( 11, -4 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-11/4" ), i.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2 3/4" ), i.str ( true ) );

    const gmp_rational j ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "9/4" ), j.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "2 1/4" ), j.str ( true ) );

    const gmp_rational k ( -18, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-9/4" ), k.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2 1/4" ), k.str ( true ) );

    const gmp_rational l ( 1, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "1/8" ), l.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "1/8" ), l.str ( true ) );

    const gmp_rational m ( 8, 1 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "8" ), m.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "8" ), m.str ( true ) );

    const gmp_rational n ( 8, 2, 1 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "10" ), n.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "10" ), n.str ( true ) );
}

void GMPTest::testIOStreamOps() {

    std::ostringstream os;
    os << gmp_rational ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "245850922/78256779" ), os.str() );

    os.str ( "" );
    os << gmp_rational ( 280.0f/375.0f );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "56/75" ), os.str() );

    std::istringstream is ( "3.14159265358979323846" );
    gmp_rational in_pi;

    is >> in_pi;

    CPPUNIT_ASSERT_EQUAL ( 245850922l, in_pi.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 78256779l, in_pi.denominator().get_si() );
}

void GMPTest::testAlgorithm() {

    const mpf_class &r (
        std::accumulate ( m_twosqrt.begin(), m_twosqrt.end(),
                          unchecked_sqrt ( 1, 1 ), std::multiplies<unchecked_sqrt > () ) );
    mp_exp_t exp;

    CPPUNIT_ASSERT_EQUAL ( std::string ( "2" ), r.get_str ( exp, 10, 4 ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           rat_vector::value_type(), std::plus<rat_vector::value_type>() )
                           .numerator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           rat_vector::value_type(), std::plus<rat_vector::value_type>() )
                           .denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           rat_vector::value_type(), std::plus<rat_vector::value_type>() )
                           .numerator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_oneseventh.begin(), m_oneseventh.end(),
                           rat_vector::value_type(),std::plus<rat_vector::value_type>() )
                           .denominator().get_si() );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
