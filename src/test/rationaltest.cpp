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
    m_accu_stein(), m_onethird(), m_accu_ul() {}

void RationalTest::setUp() {

    m_sqrt2 = Rational<uint64_t> ( std::sqrt ( 2.0 ) );

    for ( rational_type i = 1; i < 25; ++i ) {
        m_accu.push_back ( rat_vector::value_type ( 1, i ) );
        m_accu_stein.push_back ( rat_vector_stein::value_type ( 1, i ) );
    }

    m_accu_ul.reserve ( 65536u );

    for ( unsigned long ul = 1u; ul < 65536u; ++ul ) {
        m_accu_ul.push_back ( rat_vector_ul::value_type ( 1u, ul ) );
    }

    std::fill_n ( std::back_inserter ( m_onethird ), 3, rat_vector::value_type ( 1, 3 ) );
}

void RationalTest::tearDown() {}

void RationalTest::testNullRational() {
    CPPUNIT_ASSERT_EQUAL ( 0.0, static_cast<double> ( m_nullRational ) );
}

void RationalTest::testConstruct() {

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( Rational<rational_type> r ( 1, 0 ), std::runtime_error );
#endif

    Rational<rational_type, GCD_stein> a_stein ( 1, 2 );
    Rational<rational_type, GCD_stein> b_stein ( 1, -2 );
    Rational<rational_type, GCD_stein> c_stein ( -1, 2 );
    Rational<rational_type, GCD_stein> d_stein ( -1, -2 );

    CPPUNIT_ASSERT_EQUAL ( 0.5,  static_cast<double> ( a_stein ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( b_stein ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( c_stein ) );
    CPPUNIT_ASSERT_EQUAL ( 0.5,  static_cast<double> ( d_stein ) );

    CPPUNIT_ASSERT_EQUAL ( 0.5,  static_cast<double> ( Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( Rational<rational_type> ( 1, -2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( Rational<rational_type> ( -1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 0.5,  static_cast<double> ( Rational<rational_type> ( -1, -2 ) ) );

    Rational<rational_type, GCD_stein> e_stein ( 6, -8 );

    CPPUNIT_ASSERT_EQUAL ( -3, e_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4,  e_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( -3, Rational<rational_type> ( 6, -8 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4, Rational<rational_type> ( 6, -8 ).denominator() );

    Rational<rational_type, GCD_stein> f_stein ( 14, 24 );

    CPPUNIT_ASSERT_EQUAL ( 7, f_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12, f_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 7, Rational<rational_type> ( 14, 24 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12, Rational<rational_type> ( 14, 24 ).denominator() );

    Rational<rational_type, GCD_stein> g_stein ( 2, 1, 3 );

    CPPUNIT_ASSERT_EQUAL ( 7, g_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 3, g_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 7, Rational<rational_type> ( 2, 1, 3 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 3, Rational<rational_type> ( 2, 1, 3 ).denominator() );

    Rational<rational_type, GCD_stein> h_stein ( 18, 4, -5 );

    CPPUNIT_ASSERT_EQUAL ( 86, h_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5,  h_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 86, Rational<rational_type> ( 18, 4, -5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5, Rational<rational_type> ( 18, 4, -5 ).denominator() );

    Rational<rational_type, GCD_stein> i_stein ( 18, -4, 5 );

    CPPUNIT_ASSERT_EQUAL ( 86, i_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5,  i_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 86, Rational<rational_type> ( 18, -4, 5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5, Rational<rational_type> ( 18, -4, 5 ).denominator() );

    Rational<rational_type, GCD_stein> j_stein ( -18, 4, 5 );

    CPPUNIT_ASSERT_EQUAL ( -86, j_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5,   j_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( -86, Rational<rational_type> ( -18, 4, 5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5, Rational<rational_type> ( -18, 4, 5 ).denominator() );

    Rational<rational_type, GCD_stein> k_stein ( -18, 4, -5 );

    CPPUNIT_ASSERT_EQUAL ( -94, k_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5,   k_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( -94, Rational<rational_type> ( -18, 4, -5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 5, Rational<rational_type> ( -18, 4, -5 ).denominator() );
}

void RationalTest::testConstructFromDouble() {

    Rational<rational_type, GCD_euclid> p ( 19.0/51.0 );
    Rational<rational_type, GCD_stein> p_stein ( 19.0/51.0 );

    CPPUNIT_ASSERT_EQUAL ( 19, p.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 51, p.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 19, p_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 51, p_stein.denominator() );

    Rational<rational_type, GCD_euclid> q ( 516901.0/740785.0 );
    Rational<rational_type, GCD_stein> q_stein ( 516901.0/740785.0 );

    CPPUNIT_ASSERT_EQUAL ( 516901, q.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 740785, q.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 516901, q_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 740785, q_stein.denominator() );

    Rational<rational_type, GCD_euclid> r ( -0.7391304347826086 );
    Rational<rational_type, GCD_stein> r_stein ( -0.7391304347826086 );

    CPPUNIT_ASSERT_EQUAL ( -17, r.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, r.denominator() );

    CPPUNIT_ASSERT_EQUAL ( -17, r_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, r_stein.denominator() );

    Rational<rational_type, GCD_euclid> s ( 0.0 );
    Rational<rational_type, GCD_stein> s_stein ( 0.0 );

    CPPUNIT_ASSERT_EQUAL ( 0, s.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, s.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 0, s_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, s_stein.denominator() );

#ifndef __clang__
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 77227930 ), m_sqrt2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 54608393 ), m_sqrt2.denominator() );
#else
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 131836323 ), m_sqrt2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 93222358 ), m_sqrt2.denominator() );
#endif

    Rational<rational_type, GCD_euclid> pi ( M_PI );
    Rational<rational_type, GCD_stein> pi_stein ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( 245850922, pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779, pi.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 245850922, pi_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779, pi_stein.denominator() );

    CPPUNIT_ASSERT_EQUAL ( M_PI, static_cast<double> ( pi ) );
    CPPUNIT_ASSERT_EQUAL ( M_PI, static_cast<double> ( pi_stein ) );

    Rational<rational_type> t ( 1.0 );

    CPPUNIT_ASSERT_EQUAL ( 1, t.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, t.denominator() );

    Rational<rational_type> u ( 2.0 );

    CPPUNIT_ASSERT_EQUAL ( 2, u.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, u.denominator() );

    Rational<rational_type> v ( -8 );

    CPPUNIT_ASSERT_EQUAL ( -8, v.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, v.denominator() );

    Rational<uint64_t, GCD_euclid> max_pi_euclid ( 3.141592653589793238462643383279502884l );
    Rational<uint64_t, GCD_stein>  max_pi_stein ( 3.141592653589793238462643383279502884l );

    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 8717442233u ), max_pi_euclid.numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 2774848045u ), max_pi_euclid.denominator() );

    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 8717442233u ), max_pi_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint64_t> ( 2774848045u ), max_pi_stein.denominator() );
}

void RationalTest::testAssignedFromDouble() {

    Rational<rational_type> p = 19.0/51.0;

    CPPUNIT_ASSERT_EQUAL ( 19, p.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 51, p.denominator() );

    Rational<rational_type> q = 516901.0/740785.0;

    CPPUNIT_ASSERT_EQUAL ( 516901, q.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 740785, q.denominator() );

    Rational<rational_type> r = -0.7391304347826086;

    CPPUNIT_ASSERT_EQUAL ( -17, r.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, r.denominator() );

    Rational<rational_type> s = -3;

    CPPUNIT_ASSERT_EQUAL ( -3, s.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1, s.denominator() );

    Rational<rational_type> t = 1.0;

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
    Rational<rational_type, GCD_stein> pi_stein = M_PI;

    CPPUNIT_ASSERT_EQUAL ( 245850922, pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779, pi.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 245850922, pi_stein.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779, pi_stein.denominator() );

    pi = pi;

    CPPUNIT_ASSERT_EQUAL ( 245850922, pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779, pi.denominator() );
}

void RationalTest::testAddition() {

    const Rational<rational_type, GCD_euclid> a ( 17, 21 );
    const Rational<rational_type, GCD_stein> a_stein ( 17, 21 );
    const Rational<rational_type, GCD_euclid> b ( 44, 35 );
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

    const Rational<uint32_t, GCD_euclid> e ( 1, 6 );
    const Rational<uint32_t, GCD_stein> e_stein ( 1, 6 );
    const Rational<uint32_t, GCD_euclid> f ( 2, 15 );
    const Rational<uint32_t, GCD_stein> f_stein ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( 3u, ( e + f ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10u, ( e + f ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 3u, ( e_stein + f_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10u, ( e_stein + f_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 3u, ( f + e ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10u, ( f + e ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 3u, ( f_stein + e_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 10u, ( f_stein + e_stein ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 2u, ( +f ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15u, ( +f ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 2u, ( +f_stein ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 15u, ( +f_stein ).denominator() );

    Rational<rational_type> knuth_a ( 7, 66 );
    Rational<rational_type> knuth_b ( 17, 12 );

    CPPUNIT_ASSERT_EQUAL ( 67, ( knuth_a + knuth_b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( 44, ( knuth_a + knuth_b ).denominator() );
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
}

void RationalTest::testMultiplication() {

    const Rational<rational_type> a ( 2, 8 );
    const Rational<rational_type, GCD_stein> a_stein ( 2, 8 );
    const Rational<rational_type> b ( 7, 3 );
    const Rational<rational_type, GCD_stein> b_stein ( 7, 3 );

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

    CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.0f, static_cast<float> ( m_sqrt2 * m_sqrt2 ), 1.e-6 );
}

void RationalTest::testInvert() {

    CPPUNIT_ASSERT_EQUAL ( 7, Rational<rational_type> ( 161, 49 ).invert().numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, Rational<rational_type> ( 161, 49 ).invert().denominator() );

    CPPUNIT_ASSERT_EQUAL ( 7, Rational<rational_type> ( 161, 49 ).inverse().numerator() );
    CPPUNIT_ASSERT_EQUAL ( 23, Rational<rational_type> ( 161, 49 ).inverse().denominator() );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( Rational<rational_type> ().invert(), std::runtime_error );
    CPPUNIT_ASSERT_THROW ( Rational<rational_type> ().inverse(), std::runtime_error );
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
    CPPUNIT_ASSERT_THROW ( a / c, std::runtime_error );
    CPPUNIT_ASSERT_THROW ( a / ( b - d ), std::runtime_error );
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

    Rational<rational_type> d ( 542, 84 );
    Rational<rational_type, GCD_stein> d_stein ( 542, 84 );
    Rational<rational_type> e ( -65, 28 );
    Rational<rational_type, GCD_stein> e_stein ( -65, 28 );

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

    Rational<uint32_t> f ( 5.65l );
    Rational<uint32_t> g ( 1.23l );

    CPPUNIT_ASSERT_EQUAL ( static_cast<uint32_t> ( 73 ), ( f % g ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<uint32_t> ( 100 ), ( f % g ).denominator() );

    Rational<rational_type> h ( 11, 4 );

    CPPUNIT_ASSERT_EQUAL ( 2, h.mod().first );
    CPPUNIT_ASSERT_EQUAL ( 3, h.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4, h.mod().second.denominator() );

    Rational<rational_type> i ( 11, -4 );

    CPPUNIT_ASSERT_EQUAL ( -2, i.mod().first );
    CPPUNIT_ASSERT_EQUAL ( 3, i.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4, i.mod().second.denominator() );

    Rational<rational_type> j ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( 2, j.mod().first );
    CPPUNIT_ASSERT_EQUAL ( 1, j.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4, j.mod().second.denominator() );

    Rational<rational_type> k ( -18, 8 );

    CPPUNIT_ASSERT_EQUAL ( -2, k.mod().first );
    CPPUNIT_ASSERT_EQUAL ( 1, k.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4, k.mod().second.denominator() );

    Rational<rational_type> l ( 1, 8 );

    CPPUNIT_ASSERT_EQUAL ( 0, l.mod().first );
    CPPUNIT_ASSERT_EQUAL ( 1, l.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 8, l.mod().second.denominator() );

    Rational<uint32_t> m ( 18, 8 );

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

    Rational<rational_type> h ( 11, 4 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/4" ), h.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "2 3/4" ), h.str ( true ) );

    Rational<rational_type> i ( 11, -4 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-11/4" ), i.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2 3/4" ), i.str ( true ) );

    Rational<rational_type> j ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "9/4" ), j.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "2 1/4" ), j.str ( true ) );

    Rational<rational_type> k ( -18, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-9/4" ), k.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2 1/4" ), k.str ( true ) );

    Rational<rational_type> l ( 1, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "1/8" ), l.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "1/8" ), l.str ( true ) );
}

void RationalTest::testIOStreamOps() {

    std::ostringstream os;
    os << Rational<rational_type> ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "245850922/78256779" ), os.str() );

    os.str ( "" );
    os << Rational<unsigned long> ( 280.0f/375.0f );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "56/75" ), os.str() );

    std::istringstream is ( "3.14159265358979323846" );
    Rational<rational_type> in_pi;

    is >> in_pi;

    CPPUNIT_ASSERT_EQUAL ( 245850922, in_pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779, in_pi.denominator() );
}

void RationalTest::testPrecision() {

    Rational<rational_type> r ( 1, 3 );
    Rational<rational_type> s ( 2, 3 );

    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( r + r + r ) );
    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( r * 3.0 ) );
    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( 3.0 * r ) );

    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( r + s ) );
    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( s + r ) );

    Rational<rational_type> t ( -28, -963 );
    Rational<rational_type> u ( 935, 963 );

    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( t + u ) );
    CPPUNIT_ASSERT_EQUAL ( 1.0, static_cast<double> ( u + t ) );

}

void RationalTest::testAlgorithm() {

    CPPUNIT_ASSERT_DOUBLES_EQUAL ( 3.77595817775351, static_cast<double>
                                   ( std::accumulate ( m_accu.begin(), m_accu.end(),
                                           Rational<rational_type> ( 0, 1 ),
                                           std::plus<Rational<rational_type> >() ) ),
                                   m_accu.size() * std::numeric_limits<double>::epsilon() );

    CPPUNIT_ASSERT_DOUBLES_EQUAL ( 3.77595817775351, static_cast<double>
                                   ( std::accumulate ( m_accu_stein.begin(), m_accu_stein.end(),
                                           Rational<rational_type, GCD_stein> ( 0, 1 ),
                                           std::plus<Rational<rational_type, GCD_stein> >() ) ),
                                   m_accu_stein.size() * std::numeric_limits<double>::epsilon() );

    CPPUNIT_ASSERT_EQUAL ( static_cast<uint32_t> ( 3765161451u ),
                           std::accumulate ( m_accu_ul.begin(), m_accu_ul.end(),
                                   Rational<uint32_t>(),
                                   std::plus<Rational<uint32_t> >() ).numerator() );

    CPPUNIT_ASSERT_EQUAL ( 2001829888u, std::accumulate ( m_accu_ul.begin(), m_accu_ul.end(),
                           Rational<uint32_t>(), std::plus<Rational<uint32_t> >() ).denominator() );

    CPPUNIT_ASSERT_DOUBLES_EQUAL ( -3.77595817775351, static_cast<double>
                                   ( std::accumulate ( m_accu.begin(), m_accu.end(),
                                           Rational<rational_type> ( 0, 1 ),
                                           std::minus<Rational<rational_type> >() ) ),
                                   m_accu.size() * std::numeric_limits<double>::epsilon() );

    CPPUNIT_ASSERT_DOUBLES_EQUAL ( 2.08767569878681e-09,
                                   std::accumulate ( m_accu.begin(), m_accu.begin() + 12,
                                           1.0, std::multiplies<Rational<rational_type> >() ),
                                   12 * std::numeric_limits<double>::epsilon() );

    CPPUNIT_ASSERT_DOUBLES_EQUAL ( 479001600.0, std::accumulate ( m_accu.begin(),
                                   m_accu.begin() + 12, 1.0,
                                   std::divides<Rational<rational_type> >() ),
                                   12 * std::numeric_limits<double>::epsilon() );

    CPPUNIT_ASSERT_EQUAL ( 1, std::accumulate ( m_onethird.begin(),
                           m_onethird.end(), Rational<rational_type>(),
                           std::plus<Rational<rational_type> >() ).numerator() );

    CPPUNIT_ASSERT_EQUAL ( 1, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           Rational<rational_type>(),
                           std::plus<Rational<rational_type> >() ).denominator() );
}

void RationalTest::testStdMath() {

    rational_type rt;

    CPPUNIT_ASSERT_EQUAL ( std::string ( "2/3" ), std::modf ( Rational<rational_type> ( 11, 3 ) ,
                           &rt ).str() );
    CPPUNIT_ASSERT_EQUAL ( 3, rt );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), Rational<rational_type> ( 11, -3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), Rational<rational_type> ( -11, 3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), Rational<rational_type> ( 11, 3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), Rational<uint32_t> ( 11, 3 ).abs().str() );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
