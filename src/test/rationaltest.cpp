/*
 * Copyright 2015 by Heiko Schäfer <heiko@rangun.de>
 *
 * This file is part of rational.
 *
 * NetMauMau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * NetMauMau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with NetMauMau.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <numeric>
#include <sstream>

#include "rationaltest.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( RationalTest );

using namespace Commons::Math;

RationalTest::RationalTest() : CppUnit::TestFixture(), m_nullRational(), m_sqrt2(), m_accu(),
    m_accu_ul() {}

void RationalTest::setUp() {
    m_sqrt2 = Rational<rational_type> ( std::sqrt ( 2 ) );

    for ( rational_type i = 1; i < 25; ++i )
        m_accu.push_back ( rat_vector::value_type ( 1, i ) );

    m_accu_ul.reserve ( 65536u );

    for ( unsigned long ul = 1u; ul < 65536u; ++ul )
        m_accu_ul.push_back ( rat_vector_ul::value_type ( 1u, ul ) );
}

void RationalTest::tearDown() {}

void RationalTest::testNullRational() {
    CPPUNIT_ASSERT_EQUAL ( 0.0, static_cast<double> ( m_nullRational ) );
}

void RationalTest::testConstruct() {

    CPPUNIT_ASSERT_THROW ( Rational<rational_type> r ( 1, 0 ), std::runtime_error );

    CPPUNIT_ASSERT_EQUAL ( 0.5,  static_cast<double> ( Rational<rational_type> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( Rational<rational_type> ( 1, -2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( Rational<rational_type> ( -1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 0.5,  static_cast<double> ( Rational<rational_type> ( -1, -2 ) ) );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 7 ),
                           Rational<rational_type> ( 14, 24 ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 12 ),
                           Rational<rational_type> ( 14, 24 ).denominator() );
}

void RationalTest::testConstructFromDouble() {

    Rational<rational_type> p ( 19.0/51.0 );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 19 ), p.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 51 ), p.denominator() );

    Rational<rational_type> q ( 516901.0/740785.0 );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 516901 ), q.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 740785 ), q.denominator() );

    Rational<rational_type> r ( -0.7391304347826086 );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( -17 ), r.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 23 ), r.denominator() );

    Rational<rational_type> s ( 0.0 );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 0 ), s.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 1 ), s.denominator() );

#ifndef __clang__
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 77227930 ), m_sqrt2.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 54608393 ), m_sqrt2.denominator() );
#else
    CPPUNIT_ASSERT_EQUAL ( static_cast<unsigned long> ( 131836323 ), m_sqrt2.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<unsigned long> ( 93222358 ), m_sqrt2.denominator() );
#endif

    Rational<rational_type> pi ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 245850922 ), pi.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 78256779 ), pi.denominator() );

    CPPUNIT_ASSERT_EQUAL ( M_PI, static_cast<double> ( pi ) );
}

void RationalTest::testAssignedFromDouble() {

    Rational<rational_type> p = 19.0/51.0;

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 19 ), p.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 51 ), p.denominator() );

    Rational<rational_type> q = 516901.0/740785.0;

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 516901 ), q.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 740785 ), q.denominator() );

    Rational<rational_type> r = -0.7391304347826086;

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( -17 ), r.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 23 ), r.denominator() );

    Rational<rational_type> pi = M_PI;

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 245850922 ), pi.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 78256779 ), pi.denominator() );
}

void RationalTest::testAddition() {

    const Rational<rational_type> a ( 17, 21 );
    const Rational<rational_type> b ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 31 ), ( a + b ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 15 ), ( a + b ).denominator() );

    const Rational<rational_type> c ( 1, 6 );
    const Rational<rational_type> d ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 3 ), ( c + d ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 10 ), ( c + d ).denominator() );
}

void RationalTest::testSubtraction() {

    const Rational<rational_type> a ( 17, 21 );
    const Rational<rational_type> b ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( -47 ), ( a - b ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 105 ), ( a - b ).denominator() );

    const Rational<rational_type> c ( 1, 6 );
    const Rational<rational_type> d ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 1 ), ( c - d ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 30 ), ( c - d ).denominator() );
}

void RationalTest::testMultiplication() {

    const Rational<rational_type> a ( 2, 8 );
    const Rational<rational_type> b ( 7, 3 );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 7 ), ( a * b ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 12 ), ( a * b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 2.0f, static_cast<float> ( m_sqrt2 * m_sqrt2 ) );
}

void RationalTest::testInvert() {

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 7 ),
                           Rational<rational_type> ( 161, 49 ).invert().nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 23 ),
                           Rational<rational_type> ( 161, 49 ).invert().denominator() );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 7 ),
                           Rational<rational_type> ( 161, 49 ).inv().nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 23 ),
                           Rational<rational_type> ( 161, 49 ).inv().denominator() );
}

void RationalTest::testDivision() {

    const Rational<rational_type> a ( 2, 8 );
    const Rational<rational_type> b ( 7, 3 );
    const Rational<rational_type> c ( 0, 1 );
    const Rational<rational_type> d ( -7, -3 );

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 3 ), ( a / b ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 28 ), ( a / b ).denominator() );

    CPPUNIT_ASSERT_THROW ( a / c, std::runtime_error );
    CPPUNIT_ASSERT_THROW ( a / ( b - d ), std::runtime_error );
}

void RationalTest::testRelOps() {

    const Rational<rational_type> a ( 1, 4 );
    const Rational<rational_type> b ( 1, 2 );

    CPPUNIT_ASSERT ( a < b );
    CPPUNIT_ASSERT ( a <= b );

    CPPUNIT_ASSERT ( b > a );
    CPPUNIT_ASSERT ( b >= a );

    const Rational<rational_type> c ( 2, 4 );

    CPPUNIT_ASSERT ( c == b );
    CPPUNIT_ASSERT ( b == c );

    CPPUNIT_ASSERT ( a != b );
    CPPUNIT_ASSERT ( b != a );

    CPPUNIT_ASSERT ( b <= c );
    CPPUNIT_ASSERT ( c <= b );
    CPPUNIT_ASSERT ( b >= c );
    CPPUNIT_ASSERT ( c >= b );

    const Rational<rational_type> d ( 2, 4 );
    const Rational<rational_type> e ( 2, -4 );

    CPPUNIT_ASSERT ( d > e );
    CPPUNIT_ASSERT ( e < d );

    const Rational<rational_type> f ( -2, 4 );

    CPPUNIT_ASSERT ( f == e );
    CPPUNIT_ASSERT ( f >= e );
    CPPUNIT_ASSERT ( f <= e );

    CPPUNIT_ASSERT ( e == f );
    CPPUNIT_ASSERT ( e >= f );
    CPPUNIT_ASSERT ( e <= f );

    const Rational<rational_type> g ( -3, 4 );

    CPPUNIT_ASSERT ( g < d );
    CPPUNIT_ASSERT ( d > g );
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

    CPPUNIT_ASSERT ( 0.5 == Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT ( Rational<rational_type> ( 1, 2 ) == 0.5 );

    CPPUNIT_ASSERT ( 0.5 != Rational<rational_type> ( 11, 23 ) );
    CPPUNIT_ASSERT ( Rational<rational_type> ( 11, 23 ) != 0.5 );

    CPPUNIT_ASSERT ( 0.25 < Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT ( Rational<rational_type> ( 1, 2 ) > 0.25 );

    CPPUNIT_ASSERT ( 0.5 >= Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT ( Rational<rational_type> ( 1, 2 ) <= 0.5 );

    CPPUNIT_ASSERT ( 0.25 <= Rational<rational_type> ( 1, 2 ) );
    CPPUNIT_ASSERT ( Rational<rational_type> ( 1, 2 ) >= 0.25 );
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

    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 245850922 ), in_pi.nominator() );
    CPPUNIT_ASSERT_EQUAL ( static_cast<rational_type> ( 78256779 ), in_pi.denominator() );
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

    CPPUNIT_ASSERT_EQUAL ( static_cast<uint32_t> ( 104186 ), std::accumulate ( m_accu_ul.begin(),
                           m_accu_ul.end(), Rational<uint32_t>(),
                           std::plus<Rational<uint32_t> >() ).nominator() );

    CPPUNIT_ASSERT_EQUAL ( static_cast<uint32_t> ( 3502910561u ),
                           std::accumulate ( m_accu_ul.begin(), m_accu_ul.end(),
                                   Rational<uint32_t>(),
                                   std::plus<Rational<uint32_t> >() ).denominator() );

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
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
