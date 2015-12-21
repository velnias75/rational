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

#include "rationaltest.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( RationalTest );

using namespace Commons::Math;

RationalTest::RationalTest() : CppUnit::TestFixture(), m_nullRational(), m_sqrt2(), m_accu(),
    m_accu_stein(), m_onethird(), m_oneseventh(), m_accu_ul(), m_twosqrt() {}

void RationalTest::setUp() {

    m_sqrt2 = checked_sqrt ( std::sqrt ( 2.0l ) );

    for ( rational_type i = 1; i < 25; ++i ) {
        m_accu.push_back ( rat_vector::value_type ( 1, i ) );
        m_accu_stein.push_back ( rat_vector_stein::value_type ( 1, i ) );
    }

    m_accu_ul.reserve ( 47u );

    for ( rat_vector_ul::value_type::integer_type ul = 1u; ul < 47u; ++ul ) {
        m_accu_ul.push_back ( rat_vector_ul::value_type ( 1u, ul ) );
    }

    std::fill_n ( std::back_inserter ( m_onethird ), 3, rat_vector::value_type ( 1, 3 ) );
    std::fill_n ( std::back_inserter ( m_oneseventh ), 7, rat_vector::value_type ( 1, 7 ) );
    std::fill_n ( std::back_inserter ( m_twosqrt ), 2, m_sqrt2 );
}

void RationalTest::tearDown() {}

void RationalTest::testNullRational() {
    CPPUNIT_ASSERT_EQUAL ( 0.0, static_cast<double> ( m_nullRational ) );
}

void RationalTest::testConstruct() {

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( Rational<rational_type> r ( 1, 0 ), std::domain_error );
#endif

    const Rational<rational_type, GCD_stein> a_stein ( 1, 2 );
    const Rational<rational_type, GCD_stein> b_stein ( 1, -2 );
    const Rational<rational_type, GCD_stein> c_stein ( -1, 2 );
    const Rational<rational_type, GCD_stein> d_stein ( -1, -2 );

    CPPUNIT_ASSERT_EQUAL ( 0.5,  static_cast<double> ( a_stein ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( b_stein ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( c_stein ) );
    CPPUNIT_ASSERT_EQUAL ( 0.5,  static_cast<double> ( d_stein ) );

    CPPUNIT_ASSERT_EQUAL ( 0.5,  static_cast<double> ( Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( Rational<rational_type> ( 1, -2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( Rational<rational_type> ( -1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 0.5,  static_cast<double> ( Rational<rational_type> ( -1, -2 ) ) );

    const Rational<rational_type, GCD_stein> e_stein ( 6, -8 );

    CPPUNIT_ASSERT_EQUAL ( -3, e_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4,  e_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( -3, Rational<rational_type> ( 6, -8 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4, Rational<rational_type> ( 6, -8 ).denominator() );

    const Rational<rational_type, GCD_stein> f_stein ( 14, 24 );

    CPPUNIT_ASSERT_EQUAL ( 7, f_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12, f_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 7, Rational<rational_type> ( 14, 24 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12, Rational<rational_type> ( 14, 24 ).denominator() );

    const Rational<rational_type, GCD_stein> g_stein ( 2, 1, 3 );

    CPPUNIT_ASSERT_EQUAL ( 7, g_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 3, g_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 7, Rational<rational_type> ( 2, 1, 3 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 3, Rational<rational_type> ( 2, 1, 3 ).denominator() );

    const Rational<rational_type, GCD_stein> h_stein ( 18, 4, -5 );

    CPPUNIT_ASSERT_EQUAL ( 86, h_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5,  h_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 86, Rational<rational_type> ( 18, 4, -5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5, Rational<rational_type> ( 18, 4, -5 ).denominator() );

    const Rational<rational_type, GCD_stein> i_stein ( 18, -4, 5 );

    CPPUNIT_ASSERT_EQUAL ( 86, i_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5,  i_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 86, Rational<rational_type> ( 18, -4, 5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5, Rational<rational_type> ( 18, -4, 5 ).denominator() );

    const Rational<rational_type, GCD_stein> j_stein ( -18, 4, 5 );

    CPPUNIT_ASSERT_EQUAL ( -86, j_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5,   j_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( -86, Rational<rational_type> ( -18, 4, 5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5, Rational<rational_type> ( -18, 4, 5 ).denominator() );

    const Rational<rational_type, GCD_stein> k_stein ( -18, 4, -5 );

    CPPUNIT_ASSERT_EQUAL ( -94, k_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5,   k_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( -94, Rational<rational_type> ( -18, 4, -5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5, Rational<rational_type> ( -18, 4, -5 ).denominator() );
}

void RationalTest::testConstructFromDouble() {

    const Rational<rational_type, GCD_euclid> &p ( 19.0/51.0 );
    const Rational<rational_type, GCD_stein> &p_stein ( 19.0/51.0 );

    CPPUNIT_ASSERT_EQUAL ( 19, p.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 51, p.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 19, p_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 51, p_stein.denominator() );

    const Rational<rational_type, GCD_euclid> &q ( 516901.0/740785.0 );
    const Rational<rational_type, GCD_stein> &q_stein ( 516901.0/740785.0 );

    CPPUNIT_ASSERT_EQUAL ( 516901, q.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 740785, q.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 516901, q_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 740785, q_stein.denominator() );

    const Rational<rational_type, GCD_euclid> &r ( -0.7391304347826086 );
    const Rational<rational_type, GCD_stein> &r_stein ( -0.7391304347826086 );

    CPPUNIT_ASSERT_EQUAL ( -17, r.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, r.denominator() );

    CPPUNIT_ASSERT_EQUAL ( -17, r_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, r_stein.denominator() );

    const Rational<rational_type, GCD_euclid> &s ( 0.0 );
    const Rational<rational_type, GCD_stein> &s_stein ( 0.0 );

    CPPUNIT_ASSERT_EQUAL ( 0, s.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, s.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 0, s_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, s_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 6333631924u, m_sqrt2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4478554083u, m_sqrt2.denominator() );

    const Rational<rational_type, GCD_euclid> &pi ( M_PI );
    const Rational<rational_type, GCD_stein> &pi_stein ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( 245850922, pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779, pi.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 245850922, pi_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779, pi_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( M_PI, static_cast<double> ( pi ) );
    CPPUNIT_ASSERT_EQUAL ( M_PI, static_cast<double> ( pi_stein ) );

    const Rational<rational_type> &t ( 1.0 );

    CPPUNIT_ASSERT_EQUAL ( 1, t.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, t.denominator() );

    const Rational<rational_type> &u ( 2.0 );

    CPPUNIT_ASSERT_EQUAL ( 2, u.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, u.denominator() );

    const Rational<rational_type> &v ( -8 );

    CPPUNIT_ASSERT_EQUAL ( -8, v.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, v.denominator() );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( Rational<int8_t> ( 1000.0 ), std::domain_error );
#endif

    const Rational<uint64_t, GCD_euclid> max_pi_euclid ( 3.141592653589793238462643383279502884l );
    const Rational<uint64_t, GCD_stein>  max_pi_stein ( 3.141592653589793238462643383279502884l );

    CPPUNIT_ASSERT_EQUAL ( 8717442233u, max_pi_euclid.numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 2774848045u ), max_pi_euclid.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 8717442233u, max_pi_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 2774848045u ), max_pi_stein.denominator() );
}

void RationalTest::testConstructFromExpression() {

    const Rational<rational_type, GCD_euclid> &p ( "19/51" );
    const Rational<rational_type, GCD_stein> &p_stein ( "19/51" );

    CPPUNIT_ASSERT_EQUAL ( 19, p.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 51, p.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 19, p_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 51, p_stein.denominator() );

    const Rational<rational_type, GCD_euclid> &q ( "516901/740785" );
    const Rational<rational_type, GCD_stein> &q_stein ( "516901/740785" );

    CPPUNIT_ASSERT_EQUAL ( 516901, q.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 740785, q.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 516901, q_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 740785, q_stein.denominator() );

    const Rational<rational_type> &t ( "1" );

    CPPUNIT_ASSERT_EQUAL ( 1, t.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, t.denominator() );

    const Rational<rational_type> &u ( "2" );

    CPPUNIT_ASSERT_EQUAL ( 2, u.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, u.denominator() );

    const Rational<rational_type> &v ( "-8" );

    CPPUNIT_ASSERT_EQUAL ( -8, v.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, v.denominator() );

    const Rational<rational_type> &w ( "(11/2) * -8" );

    CPPUNIT_ASSERT_EQUAL ( -44, w.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, w.denominator() );

    const Rational<rational_type> &x ( "(11/2) * +(4.25+3.75)" );

    CPPUNIT_ASSERT_EQUAL ( 44, x.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, x.denominator() );

    const Rational<rational_type> &y ( "8 * -(11/2)" );

    CPPUNIT_ASSERT_EQUAL ( -44, y.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, y.denominator() );

    const Rational<rational_type> &z ( "\t8 *11.0/-2 " );

    CPPUNIT_ASSERT_EQUAL ( -44, z.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, z.denominator() );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( Rational<int8_t> ( "1000" ), std::domain_error );
#endif

    const Rational<uint64_t, GCD_euclid> max_pi_euclid ( "3.141592653589793238462643383279502884" );
    const Rational<uint64_t, GCD_stein>  max_pi_stein ( "3.141592653589793238462643383279502884" );

    CPPUNIT_ASSERT_EQUAL ( 8717442233u, max_pi_euclid.numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 2774848045u ), max_pi_euclid.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 8717442233u, max_pi_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 2774848045u ), max_pi_stein.denominator() );
}

void RationalTest::testAssignedFromDouble() {

    const Rational<rational_type> &p = 19.0/51.0;

    CPPUNIT_ASSERT_EQUAL ( 19, p.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 51, p.denominator() );

    const Rational<rational_type> &q = 516901.0/740785.0;

    CPPUNIT_ASSERT_EQUAL ( 516901, q.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 740785, q.denominator() );

    const Rational<rational_type> &r = -0.7391304347826086;

    CPPUNIT_ASSERT_EQUAL ( -17, r.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, r.denominator() );

    const Rational<rational_type> &s = -3;

    CPPUNIT_ASSERT_EQUAL ( -3, s.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, s.denominator() );

    const Rational<rational_type> &t = 1.0;

    CPPUNIT_ASSERT_EQUAL ( 1, t.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, t.denominator() );

    Rational<rational_type> u = 2.0;

    CPPUNIT_ASSERT_EQUAL ( 2, u.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, u.denominator() );

    u += 2.0;

    CPPUNIT_ASSERT_EQUAL ( 4, u.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, u.denominator() );

    u -= 2.0;

    CPPUNIT_ASSERT_EQUAL ( 2, u.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, u.denominator() );

    Rational<rational_type, GCD_euclid> pi = M_PI;
    const Rational<rational_type, GCD_stein> &pi_stein = M_PI;

    CPPUNIT_ASSERT_EQUAL ( 245850922, pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779, pi.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 245850922, pi_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779, pi_stein.denominator() );

    pi = pi;

    CPPUNIT_ASSERT_EQUAL ( 245850922, pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779, pi.denominator() );

    double v = 0.5;

    CPPUNIT_ASSERT_EQUAL ( 1.0, ( v += Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 0.5, ( v -= Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 0.25, ( v *= Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 0.5, ( v /= Rational<rational_type> ( 1, 2 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1, ( v + Rational<rational_type> ( 1, 2 ) ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( v + Rational<rational_type> ( 1, 2 ) ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 0, ( v - Rational<rational_type> ( 1, 2 ) ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( v - Rational<rational_type> ( 1, 2 ) ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( v * Rational<rational_type> ( 1, 2 ) ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4, ( v * Rational<rational_type> ( 1, 2 ) ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 2, ( v / Rational<rational_type> ( 1, 4 ) ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( v / Rational<rational_type> ( 1, 4 ) ).denominator() );
}

void RationalTest::testAddition() {

    const Rational<rational_type, GCD_euclid_fast> a ( 17, 21 );
    const Rational<rational_type, GCD_stein> a_stein ( 17, 21 );
    const Rational<rational_type, GCD_euclid_fast> b ( 44, 35 );
    const Rational<rational_type, GCD_stein> b_stein ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( 31, ( a + b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( a + b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 31, ( a_stein + b_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( a_stein + b_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 31, ( a + b_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( a + b_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 31, ( a_stein + b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( a_stein + b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 31, ( b + a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( b + a ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 31, ( b_stein + a_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( b_stein + a_stein ).denominator() );

    const Rational<rational_type, GCD_euclid> c ( 1, 6 );
    const Rational<rational_type, GCD_stein> c_stein ( 1, 6 );
    const Rational<rational_type, GCD_euclid> d ( 2, 15 );
    const Rational<rational_type, GCD_stein> d_stein ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( 3, ( c + d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10, ( c + d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 3, ( c_stein + d_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10, ( c_stein + d_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 3, ( d + c ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10, ( d + c ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 3, ( d_stein + c_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10, ( d_stein + c_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 2, ( +d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( +d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 2, ( +d_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( +d_stein ).denominator() );

    const Rational<uint32_t, GCD_euclid> e_euclid ( 1, 6 );
    const Rational<uint32_t, GCD_stein> e_stein ( 1, 6 );
    const Rational<uint32_t, GCD_euclid> f ( 2, 15 );
    const Rational<uint32_t, GCD_stein> f_stein ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( 3u, ( e_euclid + f ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10u, ( e_euclid + f ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 3u, ( e_stein + f_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10u, ( e_stein + f_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 3u, ( f + e_euclid ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10u, ( f + e_euclid ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 3u, ( f_stein + e_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10u, ( f_stein + e_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 2u, ( +f ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15u, ( +f ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 2u, ( +f_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15u, ( +f_stein ).denominator() );

    const Rational<rational_type> knuth_a ( 7, 66 );
    const Rational<rational_type> knuth_b ( 17, 12 );

    CPPUNIT_ASSERT_EQUAL ( 67, ( knuth_a + knuth_b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 44, ( knuth_a + knuth_b ).denominator() );

#ifdef __EXCEPTIONS
    const Rational<int8_t, GCD_euclid, ENABLE_OVERFLOW_CHECK> overflow ( 127, 1 );
    CPPUNIT_ASSERT_THROW ( overflow + 1.0, std::domain_error );
    const Rational<uint8_t, GCD_euclid, ENABLE_OVERFLOW_CHECK> wrap ( 255, 1 );
    CPPUNIT_ASSERT_THROW ( wrap + 1.0, std::domain_error );
#endif
}

void RationalTest::testSubtraction() {

    const Rational<rational_type, GCD_euclid> a ( 17, 21 );
    const Rational<rational_type, GCD_stein> a_stein ( 17, 21 );
    const Rational<rational_type, GCD_euclid> b ( 44, 35 );
    const Rational<rational_type, GCD_stein> b_stein ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( -47, ( a - b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 105, ( a - b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -47, ( a_stein - b_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 105, ( a_stein - b_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -47, ( a - b_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 105, ( a - b_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -47, ( a_stein - b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 105, ( a_stein - b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 0, ( a - a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( a - a ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 0, ( a_stein - a_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( a_stein - a_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 47, ( b - a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 105, ( b - a ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 47, ( b_stein - a_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 105, ( b_stein - a_stein ).denominator() );

    const Rational<rational_type, GCD_euclid> c ( 1, 6 );
    const Rational<rational_type, GCD_stein> c_stein ( 1, 6 );
    const Rational<rational_type, GCD_euclid> d ( 2, 15 );
    const Rational<rational_type, GCD_stein> d_stein ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( 1, ( c - d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 30, ( c - d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( c_stein - d_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 30, ( c_stein - d_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -1, ( d - c ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 30, ( d - c ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -1, ( d_stein - c_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 30, ( d_stein - c_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -2, ( -d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( -d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -2, ( -d_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( -d_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 2, ( d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 2, ( d_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15, ( d_stein ).denominator() );

    const Rational<int8_t> fs ( -50, 1 );

    CPPUNIT_ASSERT_EQUAL ( static_cast<int8_t> ( 50 ), ( -fs ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<int8_t> ( 1 ), ( -fs ).denominator() );

#ifdef __EXCEPTIONS
    const Rational<int8_t, GCD_euclid, ENABLE_OVERFLOW_CHECK> overflow ( -128, 1 );
    CPPUNIT_ASSERT_THROW ( overflow - 1.0, std::domain_error );
    CPPUNIT_ASSERT_THROW ( -overflow, std::domain_error );
    const Rational<uint8_t, GCD_euclid, ENABLE_OVERFLOW_CHECK> wrap ( 0, 1 );
    CPPUNIT_ASSERT_THROW ( wrap - 1.0, std::domain_error );
    CPPUNIT_ASSERT_THROW ( -wrap, std::domain_error );
#endif
}

void RationalTest::testMultiplication() {

    const Rational<rational_type> a ( 2, 8 );
    const Rational<rational_type, GCD_stein> a_stein ( 2, 8 );
    const Rational<rational_type> b ( 7, 3 );
    const Rational<rational_type, GCD_stein> b_stein ( 7, 3 );

    const Rational<rational_type> c ( -1, 1 );
    const Rational<rational_type, GCD_stein> c_stein ( -1, 1 );

    const Rational<rational_type> d ( 1, -1 );
    const Rational<rational_type, GCD_stein> d_stein ( 1, -1 );

    CPPUNIT_ASSERT_EQUAL ( 7, ( a * b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12, ( a * b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 7, ( b * a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12, ( b * a ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 7, ( a_stein * b_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12, ( a_stein * b_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 7, ( b_stein * a_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12, ( b_stein * a_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 7, ( a * b_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12, ( a * b_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 7, ( b * a_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12, ( b * a_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( c * c ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( c * c ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( c_stein * c_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( c_stein * c_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( c_stein * c ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( c_stein * c ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( c * c_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( c * c_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( d * d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( d * d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( d_stein * d_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( d_stein * d_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( d_stein * d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( d_stein * d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( d * d_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( d * d_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( c * d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( c * d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( c_stein * d_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( c_stein * d_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( c_stein * d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( c_stein * d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( c * d_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( c * d_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( d * c ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( d * c ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( d_stein * c_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( d_stein * c_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( d_stein * c ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( d_stein * c ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1, ( d * c_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( d * c_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -1, ( c * c * c ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( c * c * c ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -1, ( d * d * d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( d * d * d ).denominator() );

#ifdef __EXCEPTIONS
    const Rational<int8_t, GCD_euclid, ENABLE_OVERFLOW_CHECK> overflow ( 127, 1 );
    CPPUNIT_ASSERT_THROW ( overflow * 10.0, std::domain_error );
    const Rational<uint8_t, GCD_euclid, ENABLE_OVERFLOW_CHECK> wrap ( 255, 1 );
    CPPUNIT_ASSERT_THROW ( wrap * 2.0, std::domain_error );
    CPPUNIT_ASSERT_THROW ( m_sqrt2 * m_sqrt2, std::domain_error );
#endif
}

void RationalTest::testInvert() {

    CPPUNIT_ASSERT_EQUAL ( 7, Rational<rational_type> ( 161, 49 ).invert().numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, Rational<rational_type> ( 161, 49 ).invert().denominator() );

    CPPUNIT_ASSERT_EQUAL ( 7, Rational<rational_type> ( 161, 49 ).inverse().numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, Rational<rational_type> ( 161, 49 ).inverse().denominator() );

    CPPUNIT_ASSERT_EQUAL ( -7, Rational<rational_type> ( -161, 49 ).invert().numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, Rational<rational_type> ( -161, 49 ).invert().denominator() );

    CPPUNIT_ASSERT_EQUAL ( -7, Rational<rational_type> ( -161, 49 ).inverse().numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, Rational<rational_type> ( -161, 49 ).inverse().denominator() );

    CPPUNIT_ASSERT_EQUAL ( -7, Rational<rational_type> ( 161, -49 ).invert().numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, Rational<rational_type> ( 161, -49 ).invert().denominator() );

    CPPUNIT_ASSERT_EQUAL ( -7, Rational<rational_type> ( 161, -49 ).inverse().numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, Rational<rational_type> ( 161, -49 ).inverse().denominator() );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( Rational<rational_type> ().invert(), std::domain_error );
    CPPUNIT_ASSERT_THROW ( Rational<rational_type> ().inverse(), std::domain_error );
#endif
}

void RationalTest::testDivision() {

    const Rational<rational_type> a ( 2, 8 );
    const Rational<rational_type, GCD_stein> a_stein ( 2, 8 );
    const Rational<rational_type> b ( 7, 3 );
    const Rational<rational_type, GCD_stein> b_stein ( 7, 3 );
    const Rational<rational_type> c ( 0, 1 );
    const Rational<rational_type, GCD_stein> c_stein ( 0, 1 );
    const Rational<rational_type> d ( -7, -3 );
    const Rational<rational_type, GCD_stein> d_stein ( -7, -3 );

    CPPUNIT_ASSERT_EQUAL ( 3, ( a / b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 28, ( a / b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 28, ( b / a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 3, ( b / a ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 3, ( a_stein / b_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 28, ( a_stein / b_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 28, ( b_stein / a_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 3, ( b_stein / a_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 3, ( a / b_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 28, ( a / b_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 28, ( b / a_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 3, ( b / a_stein ).denominator() );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( a / c, std::domain_error );
    CPPUNIT_ASSERT_THROW ( a / ( b - d ), std::domain_error );
#endif
}

void RationalTest::testModulo() {

    Rational<rational_type> a ( 8, 1 );
    Rational<rational_type, GCD_stein> a_stein ( 8, 1 );

    a %= Rational<rational_type> ( 3, 1 );

    CPPUNIT_ASSERT_EQUAL ( 2, ( a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( a ).denominator() );

    a_stein %= Rational<rational_type, GCD_stein> ( 3, 1 );

    CPPUNIT_ASSERT_EQUAL ( 2, ( a_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( a_stein ).denominator() );

    a %= Rational<rational_type, GCD_stein> ( 3, 1 );

    CPPUNIT_ASSERT_EQUAL ( 2, ( a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, ( a ).denominator() );

    Rational<rational_type> c ( 41, 7 );
    c %= Rational<rational_type> ( 3, 2 );

    CPPUNIT_ASSERT_EQUAL ( 19, ( c ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 14, ( c ).denominator() );

    const Rational<rational_type> d ( 542, 84 );
    const Rational<rational_type, GCD_stein> d_stein ( 542, 84 );
    const Rational<rational_type> e ( -65, 28 );
    const Rational<rational_type, GCD_stein> e_stein ( -65, 28 );

    CPPUNIT_ASSERT_EQUAL ( -43, ( d % e ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 84, ( d % e ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 347, ( e % d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 84, ( e % d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -43, ( d_stein % e_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 84, ( d_stein % e_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 347, ( e_stein % d_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 84, ( e_stein % d_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -43, ( d % e_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 84, ( d % e_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 347, ( e % d_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 84, ( e % d_stein ).denominator() );

    const Rational<uint32_t> &f ( 5.65l );
    const Rational<uint32_t> &g ( 1.23l );

    CPPUNIT_ASSERT_EQUAL ( static_cast<uint32_t> ( 73 ), ( f % g ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint32_t> ( 100 ), ( f % g ).denominator() );

    const Rational<rational_type> h ( 11, 4 );

    CPPUNIT_ASSERT_EQUAL ( 2, h.mod().first );
    CPPUNIT_ASSERT_EQUAL ( 3, h.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4, h.mod().second.denominator() );

    const Rational<rational_type> i ( 11, -4 );

    CPPUNIT_ASSERT_EQUAL ( -2, i.mod().first );
    CPPUNIT_ASSERT_EQUAL ( -3, i.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4, i.mod().second.denominator() );

    const Rational<rational_type> j ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( 2, j.mod().first );
    CPPUNIT_ASSERT_EQUAL ( 1, j.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4, j.mod().second.denominator() );

    const Rational<rational_type> k ( -18, 8 );

    CPPUNIT_ASSERT_EQUAL ( -2, k.mod().first );
    CPPUNIT_ASSERT_EQUAL ( -1, k.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4, k.mod().second.denominator() );

    const Rational<rational_type> l ( 1, 8 );

    CPPUNIT_ASSERT_EQUAL ( 0, l.mod().first );
    CPPUNIT_ASSERT_EQUAL ( 1, l.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 8, l.mod().second.denominator() );

    const Rational<uint32_t> m ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( 2u, m.mod().first );
    CPPUNIT_ASSERT_EQUAL ( 1u, m.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4u, m.mod().second.denominator() );
}

void RationalTest::testIncDec() {

    Rational<rational_type> a ( 2, 4 );

    CPPUNIT_ASSERT_EQUAL ( 3, ( ++a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 2, ( a++ ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 5, a.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 2, a.denominator() );

    Rational<rational_type> b ( 2, 4 );

    CPPUNIT_ASSERT_EQUAL ( -1, ( --b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 2, ( b-- ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( -3, b.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 2, b.denominator() );
}

void RationalTest::testRelOps() {

    const Rational<rational_type> a ( 1, 4 );
    const Rational<rational_type, GCD_stein> a_stein ( 1, 4 );
    const Rational<rational_type> b ( 1, 2 );
    const Rational<rational_type, GCD_stein> b_stein ( 1, 2 );

    CPPUNIT_ASSERT ( a < b );
    CPPUNIT_ASSERT ( a <= b );

    CPPUNIT_ASSERT ( b > a );
    CPPUNIT_ASSERT ( b >= a );

    CPPUNIT_ASSERT ( a_stein < b_stein );
    CPPUNIT_ASSERT ( a_stein <= b_stein );

    CPPUNIT_ASSERT ( a < b_stein );
    CPPUNIT_ASSERT ( a <= b_stein );

    CPPUNIT_ASSERT ( a_stein < b );
    CPPUNIT_ASSERT ( a_stein <= b );

    CPPUNIT_ASSERT ( b_stein > a_stein );
    CPPUNIT_ASSERT ( b_stein >= a_stein );

    CPPUNIT_ASSERT ( b > a_stein );
    CPPUNIT_ASSERT ( b >= a_stein );

    CPPUNIT_ASSERT ( b_stein > a );
    CPPUNIT_ASSERT ( b_stein >= a );

    const Rational<rational_type> c ( 2, 4 );
    const Rational<rational_type, GCD_stein> c_stein ( 2, 4 );

    CPPUNIT_ASSERT ( c == b );
    CPPUNIT_ASSERT ( b == c );

    CPPUNIT_ASSERT ( c_stein == b_stein );
    CPPUNIT_ASSERT ( b_stein == c_stein );

    CPPUNIT_ASSERT ( c == b_stein );
    CPPUNIT_ASSERT ( b == c_stein );

    CPPUNIT_ASSERT ( c_stein == b );
    CPPUNIT_ASSERT ( b_stein == c );

    CPPUNIT_ASSERT ( a != b );
    CPPUNIT_ASSERT ( b != a );

    CPPUNIT_ASSERT ( a_stein != b_stein );
    CPPUNIT_ASSERT ( b_stein != a_stein );

    CPPUNIT_ASSERT ( a != b_stein );
    CPPUNIT_ASSERT ( b != a_stein );

    CPPUNIT_ASSERT ( a_stein != b );
    CPPUNIT_ASSERT ( b_stein != a );

    CPPUNIT_ASSERT ( b <= c );
    CPPUNIT_ASSERT ( c <= b );
    CPPUNIT_ASSERT ( b >= c );
    CPPUNIT_ASSERT ( c >= b );

    CPPUNIT_ASSERT ( b_stein <= c_stein );
    CPPUNIT_ASSERT ( c_stein <= b_stein );
    CPPUNIT_ASSERT ( b_stein >= c_stein );
    CPPUNIT_ASSERT ( c_stein >= b_stein );

    CPPUNIT_ASSERT ( b <= c_stein );
    CPPUNIT_ASSERT ( c <= b_stein );
    CPPUNIT_ASSERT ( b >= c_stein );
    CPPUNIT_ASSERT ( c >= b_stein );

    CPPUNIT_ASSERT ( b_stein <= c );
    CPPUNIT_ASSERT ( c_stein <= b );
    CPPUNIT_ASSERT ( b_stein >= c );
    CPPUNIT_ASSERT ( c_stein >= b );

    const Rational<rational_type> d ( 2, 4 );
    const Rational<rational_type, GCD_stein> d_stein ( 2, 4 );
    const Rational<rational_type> e ( 2, -4 );
    const Rational<rational_type, GCD_stein> e_stein ( 2, -4 );

    CPPUNIT_ASSERT ( d > e );
    CPPUNIT_ASSERT ( e < d );

    CPPUNIT_ASSERT ( d_stein > e_stein );
    CPPUNIT_ASSERT ( e_stein < d_stein );

    CPPUNIT_ASSERT ( d > e_stein );
    CPPUNIT_ASSERT ( e < d_stein );

    CPPUNIT_ASSERT ( d_stein > e );
    CPPUNIT_ASSERT ( e_stein < d );

    const Rational<rational_type> f ( -2, 4 );
    const Rational<rational_type, GCD_stein> f_stein ( -2, 4 );

    CPPUNIT_ASSERT ( f == e );
    CPPUNIT_ASSERT ( f >= e );
    CPPUNIT_ASSERT ( f <= e );

    CPPUNIT_ASSERT ( f_stein == e_stein );
    CPPUNIT_ASSERT ( f_stein >= e_stein );
    CPPUNIT_ASSERT ( f_stein <= e_stein );

    CPPUNIT_ASSERT ( f == e_stein );
    CPPUNIT_ASSERT ( f >= e_stein );
    CPPUNIT_ASSERT ( f <= e_stein );

    CPPUNIT_ASSERT ( f_stein == e );
    CPPUNIT_ASSERT ( f_stein >= e );
    CPPUNIT_ASSERT ( f_stein <= e );

    CPPUNIT_ASSERT ( e == f );
    CPPUNIT_ASSERT ( e >= f );
    CPPUNIT_ASSERT ( e <= f );

    CPPUNIT_ASSERT ( e_stein == f_stein );
    CPPUNIT_ASSERT ( e_stein >= f_stein );
    CPPUNIT_ASSERT ( e_stein <= f_stein );

    CPPUNIT_ASSERT ( e == f_stein );
    CPPUNIT_ASSERT ( e >= f_stein );
    CPPUNIT_ASSERT ( e <= f_stein );

    CPPUNIT_ASSERT ( e_stein == f );
    CPPUNIT_ASSERT ( e_stein >= f );
    CPPUNIT_ASSERT ( e_stein <= f );

    const Rational<rational_type> g ( -3, 4 );
    const Rational<rational_type, GCD_stein> g_stein ( -3, 4 );

    CPPUNIT_ASSERT ( g < d );
    CPPUNIT_ASSERT ( d > g );

    CPPUNIT_ASSERT ( g_stein < d_stein );
    CPPUNIT_ASSERT ( d_stein > g_stein );

    CPPUNIT_ASSERT ( g < d_stein );
    CPPUNIT_ASSERT ( d > g_stein );

    CPPUNIT_ASSERT ( g_stein < d );
    CPPUNIT_ASSERT ( d_stein > g );
}

void RationalTest::testGlobalOps() {

    double a = 0.5;
    a += Rational<rational_type> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 1.0, a );

    double b = a + Rational<rational_type> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 1.0, a );
    CPPUNIT_ASSERT_EQUAL ( 1.5, b );

    a -= Rational<rational_type> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.5, a );

    b = a - Rational<rational_type> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.5, a );
    CPPUNIT_ASSERT_EQUAL ( 0.0, b );

    a *= Rational<rational_type> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.25, a );

    b = a * Rational<rational_type> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.25, a );
    CPPUNIT_ASSERT_EQUAL ( 0.125, b );

    a /= Rational<rational_type> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.5, a );

    b = a / Rational<rational_type> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.5, a );
    CPPUNIT_ASSERT_EQUAL ( 1.0, b );

    double aux = 0.25;

    CPPUNIT_ASSERT_EQUAL ( 0.75, static_cast<double> ( 0.25 + Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 0.75, static_cast<double> ( Rational<rational_type> ( 1, 2 ) + 0.25 ) );
    CPPUNIT_ASSERT_EQUAL ( 0.75, aux += Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT_EQUAL ( 0.75, static_cast<double> ( Rational<rational_type> ( 1, 2 ) += 0.25 ) );

    aux = 0.25;

    CPPUNIT_ASSERT_EQUAL ( -0.25, static_cast<double> ( 0.25 - Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 0.25, static_cast<double> ( Rational<rational_type> ( 1, 2 ) - 0.25 ) );
    CPPUNIT_ASSERT_EQUAL ( -0.25, aux -= Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT_EQUAL ( 0.25, static_cast<double> ( Rational<rational_type> ( 1, 2 ) -= 0.25 ) );

    aux = 0.25;

    CPPUNIT_ASSERT_EQUAL ( 0.125, static_cast<double> ( 0.25 * Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 0.125, static_cast<double> ( Rational<rational_type> ( 1, 2 ) * 0.25 ) );
    CPPUNIT_ASSERT_EQUAL ( 0.125, aux *= Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT_EQUAL ( 0.125,
                           static_cast<double> ( Rational<rational_type> ( 1, 2 ) *= 0.25 ) );
    aux = 0.25;

    CPPUNIT_ASSERT_EQUAL ( 0.5, static_cast<double> ( 0.25 / Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 2.0, static_cast<double> ( Rational<rational_type> ( 1, 2 ) / 0.25 ) );
    CPPUNIT_ASSERT_EQUAL ( 0.5, aux /= Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT_EQUAL ( 2.0, static_cast<double> ( Rational<rational_type> ( 1, 2 ) /= 0.25 ) );

    aux = 0.25;

    CPPUNIT_ASSERT_EQUAL ( 0.25, static_cast<double> ( 0.25 % Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 0.0, static_cast<double> ( Rational<rational_type> ( 1, 2 ) % 0.25 ) );
    CPPUNIT_ASSERT_EQUAL ( 0.25, aux %= Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT_EQUAL ( 0.0, static_cast<double> ( Rational<rational_type> ( 1, 2 ) %= 0.25 ) );

    CPPUNIT_ASSERT ( 0.5 == Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT ( Rational<rational_type> ( 1, 2 ) == 0.5 );

    CPPUNIT_ASSERT ( 0.5 != Rational<rational_type> ( 11, 23 ) );
    CPPUNIT_ASSERT ( Rational<rational_type> ( 11, 23 ) != 0.5 );

    CPPUNIT_ASSERT ( 0.25 < Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT ( ! ( Rational<rational_type> ( 1, 2 ) < 0.25 ) );

    CPPUNIT_ASSERT ( ! ( 0.25 > Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT ( Rational<rational_type> ( 1, 2 ) > 0.25 );

    CPPUNIT_ASSERT ( 0.5 >= Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT ( Rational<rational_type> ( 1, 2 ) <= 0.5 );

    CPPUNIT_ASSERT ( 0.25 <= Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT ( Rational<rational_type> ( 1, 2 ) >= 0.25 );

    CPPUNIT_ASSERT ( ( !Rational<rational_type> ( 1, 2 ) ) == false );
    CPPUNIT_ASSERT ( ( !Rational<rational_type> ( 0, 2 ) ) == true );
    CPPUNIT_ASSERT ( ( !Rational<rational_type> ( 0, -2 ) ) == true );
}

void RationalTest::testString() {

    const Rational<rational_type> h ( 11, 4 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/4" ), h.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "2 3/4" ), h.str ( true ) );

    const Rational<rational_type> i ( 11, -4 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-11/4" ), i.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2 3/4" ), i.str ( true ) );

    const Rational<rational_type> j ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "9/4" ), j.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "2 1/4" ), j.str ( true ) );

    const Rational<rational_type> k ( -18, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-9/4" ), k.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2 1/4" ), k.str ( true ) );

    const Rational<rational_type> l ( 1, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "1/8" ), l.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "1/8" ), l.str ( true ) );

    const Rational<rational_type> m ( 8, 1 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "8" ), m.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "8" ), m.str ( true ) );

    const Rational<rational_type> n ( 8, 2, 1 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "10" ), n.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "10" ), n.str ( true ) );
}

void RationalTest::testIOStreamOps() {

    std::istringstream real_in ( "0.333333333" );
    Rational<rational_type, GCD_euclid_fast, ENABLE_OVERFLOW_CHECK> real_rat;

    real_in >> real_rat;

    CPPUNIT_ASSERT_EQUAL ( 333333333, real_rat.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1000000000, real_rat.denominator() );

    std::istringstream expr_in ( "1/8 * 8.897 - 3.6 *\n1" );
    Rational<rational_type> expr_rat;

    expr_in >> expr_rat;

    CPPUNIT_ASSERT_EQUAL ( -19903, expr_rat.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 8000, expr_rat.denominator() );

    std::ostringstream os;
    os << Rational<rational_type> ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "245850922/78256779" ), os.str() );

    os.str ( "" );
    os << Rational<unsigned long> ( 280.0f/375.0f );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "56/75" ), os.str() );

    Rational<uint64_t> in_pi;

    std::istringstream is ( "3.14159265358979323846" );
    is >> in_pi;

    CPPUNIT_ASSERT_EQUAL ( 8717442233u, in_pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 2774848045u ), in_pi.denominator() );

    os.str ( "" );
    os << in_pi;

    std::istringstream is2 ( os.str() );
    is2 >> in_pi;

    CPPUNIT_ASSERT_EQUAL ( 8717442233u, in_pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 2774848045u ), in_pi.denominator() );

#ifdef __EXCEPTIONS
    const Rational<uint64_t> out_mixed ( 2, 3, 4 );

    os.str ( "" );
    os << out_mixed.str ( true );

    std::istringstream is3 ( os.str() );

    CPPUNIT_ASSERT_THROW ( is3 >> in_pi, std::runtime_error );
#endif
}

void RationalTest::testPrecision() {

    const Rational<rational_type> r ( 1, 3 );
    const Rational<rational_type> s ( 2, 3 );

    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( r + r + r ) );
    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( r * 3.0 ) );
    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( 3.0 * r ) );

    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( r + s ) );
    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( s + r ) );

    const Rational<rational_type> t ( -28, -963 );
    const Rational<rational_type> u ( 935, 963 );

    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( t + u ) );
    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( u + t ) );

}

#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic push
void RationalTest::testAlgorithm() {

    CPPUNIT_ASSERT_DOUBLES_EQUAL ( 3.77595817775351, static_cast<double>
                                   ( std::accumulate ( m_accu.begin(), m_accu.end(),
                                           rat_vector::value_type ( 0, 1 ),
                                           std::plus<rat_vector::value_type >() ) ),
                                   m_accu.size() * std::numeric_limits<double>::epsilon() );

    CPPUNIT_ASSERT_DOUBLES_EQUAL ( 3.77595817775351, static_cast<double>
                                   ( std::accumulate ( m_accu_stein.begin(), m_accu_stein.end(),
                                           Rational<rat_vector_stein::value_type::integer_type,
                                           GCD_stein> ( 0, 1 ), std::plus<
                                           Rational<rat_vector_stein::value_type::integer_type,
                                           GCD_stein> >() ) ),
                                   m_accu_stein.size() * std::numeric_limits<double>::epsilon() );

    CPPUNIT_ASSERT_EQUAL ( 5943339269060627227u, std::accumulate ( m_accu_ul.begin(),
                           m_accu_ul.end(), rat_vector_ul::value_type(),
                           std::plus<rat_vector_ul::value_type>() ).numerator() );

    CPPUNIT_ASSERT_EQUAL ( 1345655451257488800u, std::accumulate ( m_accu_ul.begin(),
                           m_accu_ul.end(), rat_vector_ul::value_type(),
                           std::plus<rat_vector_ul::value_type>() ).denominator() );

    CPPUNIT_ASSERT_DOUBLES_EQUAL ( -3.77595817775351, static_cast<double>
                                   ( std::accumulate ( m_accu.begin(), m_accu.end(),
                                           rat_vector::value_type ( 0, 1 ),
                                           std::minus<rat_vector::value_type>() ) ),
                                   m_accu.size() * std::numeric_limits<double>::epsilon() );

    CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.08767569878681e-09,
                                   std::accumulate ( m_accu.begin(), m_accu.begin() + 12,
                                           1.0, std::multiplies<rat_vector::value_type>() ),
                                   12.0 * std::numeric_limits<double>::epsilon() );

    CPPUNIT_ASSERT_DOUBLES_EQUAL ( 479001600.0, std::accumulate ( m_accu.begin(),
                                   m_accu.begin() + 12, 1.0,
                                   std::divides<rat_vector::value_type>() ),
                                   12.0 * std::numeric_limits<double>::epsilon() );

    CPPUNIT_ASSERT_EQUAL ( 1, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           rat_vector::value_type(), std::plus<rat_vector::value_type>() )
                           .numerator() );

    CPPUNIT_ASSERT_EQUAL ( 1, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           rat_vector::value_type(), std::plus<rat_vector::value_type>() ).
                           denominator() );

    CPPUNIT_ASSERT_EQUAL ( 1.0l, static_cast<long double> ( std::accumulate ( m_onethird.begin(),
                           m_onethird.end(), rat_vector::value_type(),
                           std::plus<rat_vector::value_type>() ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1.0l, static_cast<long double> ( std::accumulate ( m_oneseventh.begin(),
                           m_oneseventh.end(), rat_vector::value_type(),
                           std::plus<rat_vector::value_type>() ) ) );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( std::accumulate ( m_twosqrt.begin(), m_twosqrt.end(),
                           checked_sqrt ( 1, 1 ), std::multiplies<checked_sqrt > () ),
                           std::domain_error );
#endif

    const Rational<rational_type> a ( 77, 88 );
    const Rational<rational_type> b ( 88, 77 );

    CPPUNIT_ASSERT_EQUAL ( a, std::min ( a, b ) );
    CPPUNIT_ASSERT_EQUAL ( a, std::min ( b, a ) );

    CPPUNIT_ASSERT_EQUAL ( b, std::max ( a, b ) );
    CPPUNIT_ASSERT_EQUAL ( b, std::max ( b, a ) );

    const rational_type cf_pi[5] = { 3, 7, 15, 1, 292 };

    CPPUNIT_ASSERT_EQUAL ( 103993, cf ( cf_pi, cf_pi + 5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 33102, cf ( cf_pi, cf_pi + 5 ).denominator() );

    std::vector<rational_type> o_pi;
    seq ( cf ( cf_pi, cf_pi + 5 ), std::back_inserter ( o_pi ) );

    CPPUNIT_ASSERT_EQUAL ( static_cast<std::vector<rational_type>::size_type> ( 5u ), o_pi.size() );
    CPPUNIT_ASSERT ( std::equal ( o_pi.begin(), o_pi.end(), cf_pi ) );

    rational_type ccf[] = { 0, 3 };
    std::vector<rational_type> ocf;

    seq ( Rational<rational_type> ( 1, 3 ), std::back_inserter ( ocf ) );

    CPPUNIT_ASSERT_EQUAL ( static_cast<std::vector<rational_type>::size_type> ( 2u ), ocf.size() );
    CPPUNIT_ASSERT ( std::equal ( ocf.begin(), ocf.end(), ccf ) );

    const Rational<rational_type> c ( 88, -77 );
    const rational_type ancf[] = { -1, -7 };
    std::vector<rational_type> negcf;

    seq ( c, std::back_inserter ( negcf ) );

    CPPUNIT_ASSERT ( std::equal ( negcf.begin(), negcf.end(), ancf ) );
}
#pragma GCC diagnostic pop

void RationalTest::testStdMath() {

    rational_type rt;

    CPPUNIT_ASSERT_EQUAL ( std::string ( "2/3" ), std::modf ( Rational<rational_type> ( 11, 3 ) ,
                           &rt ).str() );
    CPPUNIT_ASSERT_EQUAL ( 3, rt );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), Rational<rational_type> ( 11, -3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), Rational<rational_type> ( -11, 3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), Rational<rational_type> ( 11, 3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), Rational<uint32_t> ( 11, 3 ).abs().str() );

    const Rational<rational_type> &a ( Rational<rational_type>::rf_info ( 142857 ) );

    CPPUNIT_ASSERT_EQUAL ( 1, a.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 7, a.denominator() );

    const Rational<rational_type> &b ( Rational<rational_type>::rf_info ( 34 ) );

    CPPUNIT_ASSERT_EQUAL ( 34, b.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 99, b.denominator() );

    const Rational<rational_type> &c ( Rational<rational_type>::rf_info ( 123456789 ) );

    CPPUNIT_ASSERT_EQUAL ( 13717421, c.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 111111111, c.denominator() );

    const Rational<rational_type> &d ( Rational<rational_type>::rf_info ( 12, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( 4, d.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 333, d.denominator() );

    const Rational<rational_type> &ex ( Rational<rational_type>::rf_info ( 6, 0, 1111 ) );

    CPPUNIT_ASSERT_EQUAL ( 667, ex.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 6000, ex.denominator() );

    const Rational<unsigned long> &f ( Rational<unsigned long>::rf_info ( 1, 2, 3, 4 ) );

    CPPUNIT_ASSERT_EQUAL ( 1499ul, f.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 49950000ul, f.denominator() );

    const Rational<unsigned long> &g ( Rational<unsigned long>::rf_info ( 6, 0, 0, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( 1ul, g.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15ul, g.denominator() );

    const Rational<unsigned long> &h ( Rational<unsigned long>::rf_info ( 6, 0, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( 1ul, h.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 6ul, h.denominator() );

    const Rational<unsigned long> &i ( Rational<unsigned long>::rf_info ( 1, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( 1ul, i.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 99ul, i.denominator() );

    const Rational<unsigned long> &j ( Rational<unsigned long>::rf_info ( 1 ) );

    CPPUNIT_ASSERT_EQUAL ( 1ul, j.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 9ul, j.denominator() );

    Rational<unsigned long>::rf_info dc;

    const Rational<unsigned long> k ( 7ul, 13ul );

    const unsigned long k_digits[] = { 5, 3, 8, 4, 6, 1 };

    CPPUNIT_ASSERT_EQUAL ( 0ul, k.decompose ( dc ) );
    CPPUNIT_ASSERT_EQUAL ( 7ul, Rational<unsigned long> ( dc ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 13ul, Rational<unsigned long> ( dc ).denominator() );

    CPPUNIT_ASSERT ( std::equal ( dc.reptend_digits.begin(), dc.reptend_digits.begin(),
                                  k_digits ) );

    const Rational<unsigned long> l ( 88ul, 100ul );

    CPPUNIT_ASSERT_EQUAL ( 0ul, l.decompose ( dc ) );
    CPPUNIT_ASSERT_EQUAL ( 22ul, Rational<unsigned long> ( dc ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 25ul, Rational<unsigned long> ( dc ).denominator() );

    const Rational<unsigned long> m ( 8ul, 3ul );

    CPPUNIT_ASSERT_EQUAL ( 2ul, m.decompose ( dc ) );
    CPPUNIT_ASSERT_EQUAL ( 2ul, Rational<unsigned long> ( dc ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 3ul, Rational<unsigned long> ( dc ).denominator() );

    const Rational<unsigned long> n ( "(70/2) - (1741832/249975)" );

    CPPUNIT_ASSERT_EQUAL ( 28ul, n.decompose ( dc ) );
    CPPUNIT_ASSERT_EQUAL ( 3ul, dc.pre );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 1 ), dc.pre_leading_zeros );
    CPPUNIT_ASSERT_EQUAL ( 1975ul, dc.reptend );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 0 ), dc.leading_zeros );

    Rational<long>::rf_info sdc;

    const Rational<long> o ( -3, 1, 3 );

    CPPUNIT_ASSERT_EQUAL ( -2l, o.decompose ( sdc ) );
    CPPUNIT_ASSERT ( sdc.pre_digits.empty() );
    CPPUNIT_ASSERT_EQUAL ( -6l, sdc.reptend_digits.front() );

    const Rational<long> p ( 13, -30 );

    CPPUNIT_ASSERT_EQUAL ( 0l, p.decompose ( sdc ) );
    CPPUNIT_ASSERT_EQUAL ( -4l, sdc.pre_digits.front() );
    CPPUNIT_ASSERT_EQUAL ( -3l, sdc.reptend_digits.front() );

    const Rational<long> q ( -2, 5 );

    CPPUNIT_ASSERT_EQUAL ( 0l, q.decompose ( sdc ) );
    CPPUNIT_ASSERT_EQUAL ( -4l, sdc.pre_digits.front() );
    CPPUNIT_ASSERT ( sdc.reptend_digits.empty() );

    const Rational<long> r ( 8, -2, 5 );

    CPPUNIT_ASSERT_EQUAL ( 7l, r.decompose ( sdc ) );
    CPPUNIT_ASSERT_EQUAL ( 6l, sdc.pre_digits.front() );
    CPPUNIT_ASSERT ( sdc.reptend_digits.empty() );

    const Rational<unsigned long> s ( 3, 4 );

    CPPUNIT_ASSERT_EQUAL ( 243ul, s.pow ( 5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1024ul, s.pow ( 5 ).denominator() );

#ifdef __EXCEPTIONS
    const Rational<long> t ( 3, 4 );

    CPPUNIT_ASSERT_THROW ( t.pow ( 0 ), std::domain_error );
    CPPUNIT_ASSERT_THROW ( t.pow ( -8 ), std::domain_error );
#endif
}

void RationalTest::testRatRat() {

    const Rational<rational_type> a ( 77, 88 );
    const Rational<rational_type> b ( 88, 77 );
    const Rational<rational_type> c ( a, b );

    CPPUNIT_ASSERT_EQUAL ( 49, c.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 64, c.denominator() );

    const Rational<uint32_t, GCD_euclid,      ENABLE_OVERFLOW_CHECK> d ( 7, 8 );
    const Rational<uint32_t, GCD_euclid_fast, NO_OPERATOR_CHECK>     e ( 8, 7 );
    const Rational<uint32_t, GCD_stein,       ENABLE_OVERFLOW_CHECK> f ( d, e );

    CPPUNIT_ASSERT_EQUAL ( 49u, f.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 64u, f.denominator() );

    const Rational<rational_type> g ( Rational<rational_type> ( 88 ), a );

    CPPUNIT_ASSERT_EQUAL ( 704, g.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 7, g.denominator() );
}

void RationalTest::testGoldenRatio() {

    Rational<uint64_t, GCD_null> phi ( 1u, 1u );

    for ( std::size_t i = 0u; i < 91u; ++i ) ( ++phi ).invert();

    CPPUNIT_ASSERT_EQUAL ( 12200160415121876738u, phi.inverse().numerator() );
    CPPUNIT_ASSERT_EQUAL ( 7540113804746346429u, phi.inverse().denominator() );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
