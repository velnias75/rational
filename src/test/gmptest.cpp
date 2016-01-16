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

#include <numeric>

#include "gmptest.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( GMPTest );

using namespace Commons::Math;

GMPTest::GMPTest() : CppUnit::TestFixture(), m_sqrt2(), m_twosqrt(), m_onethird(),
    m_oneseventh() {}

void GMPTest::setUp() {

    using namespace std;

    m_sqrt2 = unchecked_sqrt ( mpf_class ( sqrt ( mpf_class ( 2.0 ) ) ) );

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
        std::string ( "4478554083" ), m_sqrt2.numerator().get_str() );

    CPPUNIT_ASSERT_EQUAL (
        std::string ( "3166815962" ), m_sqrt2.denominator().get_str() );

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

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-1331499020266063" ), r.numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "1801439850948203" ), r.denominator().get_str() );

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
    CPPUNIT_ASSERT_EQUAL ( -3l, i.mod().second.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4l, i.mod().second.denominator().get_si() );

    gmp_rational j ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( 2l, j.mod().first.get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, j.mod().second.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4l, j.mod().second.denominator().get_si() );

    gmp_rational k ( -18, 8 );

    CPPUNIT_ASSERT_EQUAL ( -2l, k.mod().first.get_si() );
    CPPUNIT_ASSERT_EQUAL ( -1l, k.mod().second.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4l, k.mod().second.denominator().get_si() );

    gmp_rational l ( 1, 8 );

    CPPUNIT_ASSERT_EQUAL ( 0l, l.mod().first.get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, l.mod().second.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 8l, l.mod().second.denominator().get_si() );

    const gmp_rational n ( 2, 1 );

    CPPUNIT_ASSERT_EQUAL ( 2l, n.mod().first.get_si() );
    CPPUNIT_ASSERT_EQUAL ( 0l, n.mod().second.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, n.mod().second.denominator().get_si() );
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

    std::istringstream real_in ( "0.33333333333333333" );
    gmp_rational real_rat;

    real_in >> real_rat;

#ifdef HAVE_MPREAL_H
    CPPUNIT_ASSERT_EQUAL ( std::string ( "1" ), real_rat.numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "3" ), real_rat.denominator().get_str() );
#else
    CPPUNIT_ASSERT_EQUAL ( std::string ( "33357584220089605" ), real_rat.numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "100072752660268816" ), real_rat.denominator().get_str() );
#endif

    std::ostringstream os;
    os << gmp_rational ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "245850922/78256779" ), os.str() );

    os.str ( "" );
    os << gmp_rational ( 280.0f/375.0f );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "56/75" ), os.str() );

    gmp_rational in_pi;
    std::istringstream is ( "(3 + 0.14159265358979323846)" );

    is >> in_pi;

#ifdef HAVE_MPREAL_H
    CPPUNIT_ASSERT_EQUAL ( std::string ( "657408909" ), in_pi.numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "209259755" ), in_pi.denominator().get_str() );
#else
    CPPUNIT_ASSERT_EQUAL ( std::string ( "21053343141" ), in_pi.numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "6701487259" ), in_pi.denominator().get_str() );
#endif

}

#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic push
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

    const mpz_class cf_pi[] = {  3,  7, 15,  1, 292,  1, 1,   1,  2,  1,  3, 1, 14, 2, 1,  1, 2, 2,
                                 2,  2,  1, 84,   2,  1, 1,  15,  3, 13,  1, 4,  2, 6, 6, 99, 1, 2,
                                 2,  6,  3,  5,   1,  1, 6,   8,  1,  7,  1, 2,  3, 7, 1,  2, 1, 1,
                                 12, 1,  1,  1,   3,  1, 1,   8,  1,  1,  2, 1,  6, 1, 1,  5, 2, 2,
                                 3,  1,  2,  4,   4, 16, 1, 161, 45,  1, 22, 1,  2, 2, 1,  4, 1, 2,
                                 24, 1,  2,  1,   3,  1, 3
                              };

    CPPUNIT_ASSERT_EQUAL ( std::string ( "227159758552934520439668309319746303422708645581861" ),
                           cf ( cf_pi, cf_pi + 97 ).numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "72307196890521956737416455481060519150048966236850" ),
                           cf ( cf_pi, cf_pi + 97 ).denominator().get_str() );

    std::vector<gmp_rational::integer_type> o_pi;
    seq ( cf ( cf_pi, cf_pi + 97 ), std::back_inserter ( o_pi ) );

    CPPUNIT_ASSERT_EQUAL ( static_cast<std::vector<gmp_rational::integer_type>::size_type> ( 97u ),
                           o_pi.size() );
    CPPUNIT_ASSERT ( std::equal ( o_pi.begin(), o_pi.end(), cf_pi ) );
}
#pragma GCC diagnostic pop

void GMPTest::testStdMath() {

    rational_type rt;

    CPPUNIT_ASSERT_EQUAL ( std::string ( "2/3" ), std::modf ( gmp_rational ( 11, 3 ) ,
                           &rt ).str() );
    CPPUNIT_ASSERT_EQUAL ( rational_type ( 3 ), rt );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), gmp_rational ( 11, -3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), gmp_rational ( -11, 3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), gmp_rational ( 11, 3 ).abs().str() );

    const gmp_rational &a ( gmp_rational::rf_info ( 142857 ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, a.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 7l, a.denominator().get_si() );

    const gmp_rational &b ( gmp_rational::rf_info ( 34 ) );

    CPPUNIT_ASSERT_EQUAL ( 34l, b.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 99l, b.denominator().get_si() );

    const gmp_rational &c ( gmp_rational::rf_info ( 123456789 ) );

    CPPUNIT_ASSERT_EQUAL ( 13717421l, c.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 111111111l, c.denominator().get_si() );

    const gmp_rational &d ( gmp_rational::rf_info ( 12, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( 4l, d.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 333l, d.denominator().get_si() );

    const gmp_rational &ex ( gmp_rational::rf_info ( 6, 0, 1111 ) );

    CPPUNIT_ASSERT_EQUAL ( 667l, ex.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 6000l, ex.denominator().get_si() );

    const gmp_rational &f ( gmp_rational::rf_info ( 1, 2, 3, 4 ) );

    CPPUNIT_ASSERT_EQUAL ( 1499l, f.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 49950000l, f.denominator().get_si() );

    const gmp_rational &g ( gmp_rational::rf_info ( 6, 0, 0, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, g.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 15l, g.denominator().get_si() );

    const gmp_rational &h ( gmp_rational::rf_info ( 6, 0, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, h.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 6l, h.denominator().get_si() );

    const gmp_rational &i ( gmp_rational::rf_info ( 1, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, i.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 99l, i.denominator().get_si() );

    const gmp_rational &j ( gmp_rational::rf_info ( 1 ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, j.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 9l, j.denominator().get_si() );

    gmp_rational::rf_info dc;
    std::vector<gmp_rational::rf_info::digit_type> pre, rep;

    const gmp_rational k ( 7, 13 );

    CPPUNIT_ASSERT_EQUAL ( 0l, k.decompose ( dc, pre, rep ).get_si() );
    CPPUNIT_ASSERT_EQUAL ( 7l, gmp_rational ( dc ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 13l, gmp_rational ( dc ).denominator().get_si() );

    const gmp_rational l ( 88, 100 );

    CPPUNIT_ASSERT_EQUAL ( 0l, l.decompose ( dc, pre, rep ).get_si() );
    CPPUNIT_ASSERT_EQUAL ( 22l, gmp_rational ( dc ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 25l, gmp_rational ( dc ).denominator().get_si() );

    const gmp_rational m ( 8, 3 );

    CPPUNIT_ASSERT_EQUAL ( 2l, m.decompose ( dc, pre, rep ).get_si() );
    CPPUNIT_ASSERT_EQUAL ( 2l, gmp_rational ( dc ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 3l, gmp_rational ( dc ).denominator().get_si() );

    const gmp_rational n ( 8, 1, 53 );

    const gmp_rational::integer_type n_digits[] = { 0, 1, 8, 8, 6, 7, 9, 2, 4, 5, 2, 8, 3 };

    CPPUNIT_ASSERT_EQUAL ( 8l, n.decompose ( dc, pre, rep ).get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, gmp_rational ( dc ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 53l, gmp_rational ( dc ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "188679245283" ), dc.reptend.get_str() );
    CPPUNIT_ASSERT ( std::equal ( rep.begin(), rep.end(), n_digits ) );

    const gmp_rational o ( 1, 31 );

    CPPUNIT_ASSERT_EQUAL ( 0l, o.decompose ( dc, pre, rep ).get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, gmp_rational ( dc ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 31l, gmp_rational ( dc ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "32258064516129" ), dc.reptend.get_str() );

    const gmp_rational p ( "(70/2) - (1741832/249975)" );

    CPPUNIT_ASSERT_EQUAL ( 28l, p.decompose ( dc, pre, rep ).get_si() );
    CPPUNIT_ASSERT_EQUAL ( 3l, dc.pre.get_si() );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 1 ), dc.pre_leading_zeros );
    CPPUNIT_ASSERT_EQUAL ( 1975l, dc.reptend.get_si() );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 0 ), dc.leading_zeros );

    const gmp_rational q ( "123.32 / (12453/370)" );

    CPPUNIT_ASSERT_EQUAL ( 228142l, q.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 62265l, q.denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 3l, q.decompose ( dc, pre, rep ).get_si() );
    CPPUNIT_ASSERT_EQUAL ( 6l, dc.pre.get_si() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "64048823576648197221553039428250220830322010760459327" \
                                         "06978238175540030514735405123263470649642656388018951" \
                                         "25672528707941861398859712519071709628202039669156026" \
                                         "66024251184453545330442463663374287320324419818517626" \
                                         "27479322251666265156990283465831526539789608929575202" \
                                         "76238657351642174576407291415723118927166144704087368" \
                                         "50558098450172649160844776359110254557134826949329478" \
                                         "84044005460531598811531357905725527985224443909098209" \
                                         "26684333092427527503412832249257207098691078454990765" \
                                         "27744318638079177708182767204689633020155785754436681" \
                                         "92403436922829840199148799486067614229502931020637597" \
                                         "36609652292620252148076768650124467999678792258893439" \
                                         "33188789849835381032682887657592547980406327792499799" \
                                         "24516180839958242993656147113145426804785995342487753" \
                                         "95487031237452822613024973901871035091945715891752991" \
                                         "24708905484622179394523408014133140608688669396932466" \
                                         "07243234561952943065927888862121577130008833212880430" \
                                         "41837308279129527021601220589416204930538825985706255" \
                                         "52075805026901148317674455954388500762868385128081586" \
                                         "76624106640970047378141813217698546534971492812976792" \
                                         "74070505099172890066650606279611338633261061591584357" \
                                         "18300811049546294065686983056291656628924757086645788" \
                                         "16349474022323938006905966433791054364410182285393077" \
                                         "97317915361760218421263952461254316229021119408977756" \
                                         "36392837067373323697101100136513289970288283947643138" \
                                         "19963061109772745523167108327310688187585320806231430" \
                                         "17746727696137476913193607965951979442704569180117240" \
                                         "82550389464386091704810085923070746004978719987151690" \
                                         "35573757327551593993415241307315506303701919216253111" \
                                         "69999196980647233598329719746245884525817072191439813" \
                                         "69951015819481249498112904520998956074841403677828635" \
                                         "67011964988356219384887175780936320565325624347546775" \
                                         "87729864289729382478117722637115554484863085200353328" \
                                         "515217216734923311651810808" ), dc.reptend.get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 0 ), dc.pre_leading_zeros );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 1776 ), rep.size() );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 0 ), dc.leading_zeros );

    const gmp_rational s ( 3, 4 );

    CPPUNIT_ASSERT_EQUAL ( 81l, s.pow ( 4 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 256l, s.pow ( 4 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 243l, s.pow ( 5 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1024l, s.pow ( 5 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "485192780976896426811558553967593360" \
                                         "72749841943521979872827" ),
                           s.pow ( 123 ).numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "113078212145816597093331040047546785" \
                                         "012958969400039613319782796882727665" \
                                         "664" ), s.pow ( 123 ).denominator().get_str() );

#ifdef __EXCEPTIONS
    const gmp_rational t ( 3, 4 );

    CPPUNIT_ASSERT_THROW ( t.pow ( 0 ), std::domain_error );
    CPPUNIT_ASSERT_THROW ( t.pow ( -8 ), std::domain_error );
#endif

    const gmp_rational u ( 2, 1 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "4946041176255201878775086487573351061418968498177" ),
                           u.sqrt().numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "3497379255757941172020851852070562919437964212608" ),
                           u.sqrt().denominator().get_str() );

    const gmp_rational v ( 10, 17 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "1983567417147843927170789761" ),
                           v.sqrt().numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "2586255495350365951590026592" ),
                           v.sqrt().denominator().get_str() );

    const gmp_rational w ( 9, 1 );

    CPPUNIT_ASSERT_EQUAL ( 3l, w.sqrt().numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, w.sqrt().denominator().get_si() );

    const gmp_rational x ( mpz_class ( "785791622400625" ), 1 );

    CPPUNIT_ASSERT_EQUAL ( 28031975l, x.sqrt().numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, x.sqrt().denominator().get_si() );

    const gmp_rational y ( 256, 81 );

    CPPUNIT_ASSERT_EQUAL ( 16l, y.sqrt().numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 9l, y.sqrt().denominator().get_si() );

#ifdef __EXCEPTIONS
    const gmp_rational z ( -256, 81 );
    CPPUNIT_ASSERT_THROW ( z.sqrt(), std::domain_error );
#endif
}

void GMPTest::testGoldenRatio() {

    Rational<gmp_rational::integer_type, GCD_null> phi ( gmp_rational::one_, gmp_rational::one_ );

    for ( std::size_t i = 0u; i < 1024u; ++i ) ( ++phi ).invert();

    CPPUNIT_ASSERT_EQUAL ( std::string ( "1179869281805523255014757888412586560808902854456091" \
                                         "3468519228968187430794620907976123201977895385245239" \
                                         "7050828306569046301783141598663704952115390234610526" \
                                         "8281123032179655593090772272438413164852733945840731" \
                                         "7543768" ), phi.inverse().numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "7291993184377412737043195648396979558721167948342308" \
                                         "6377162058185874001489121865798744093687543548489948" \
                                         "3181625031189341064810479244078947534047137736685242" \
                                         "0526027975140687031196633477605718294523235826853392" \
                                         "138525" ), phi.inverse().denominator().get_str() );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
