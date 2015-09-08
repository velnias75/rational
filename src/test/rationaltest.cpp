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

#include <sstream>

#include "rationaltest.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( RationalTest );

using namespace Commons::Math;

void RationalTest::setUp() {
    m_sqrt2 = Rational<long long> ( std::sqrt ( 2 ) );
}

void RationalTest::tearDown() {}

void RationalTest::testNullRational() {
    CPPUNIT_ASSERT_EQUAL ( 0.0, static_cast<double> ( m_nullRational ) );
}

void RationalTest::testConstruct() {

    CPPUNIT_ASSERT_EQUAL ( 0.5, static_cast<double> ( Rational<long long> ( 1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( Rational<long long> ( 1, -2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, static_cast<double> ( Rational<long long> ( -1, 2 ) ) );
    CPPUNIT_ASSERT_EQUAL ( 0.5, static_cast<double> ( Rational<long long> ( -1, -2 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 7ll,  Rational<long long> ( 14, 24 ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( 12ll, Rational<long long> ( 14, 24 ).denominator() );
}

void RationalTest::testConstructFromDouble() {

    Rational<long long> p ( 19.0/51.0 );

    CPPUNIT_ASSERT_EQUAL ( 19ll, p.nominator() );
    CPPUNIT_ASSERT_EQUAL ( 51ll, p.denominator() );

    Rational<long long> q ( 516901.0/740785.0 );

    CPPUNIT_ASSERT_EQUAL ( 516901ll, q.nominator() );
    CPPUNIT_ASSERT_EQUAL ( 740785ll, q.denominator() );

    Rational<long long> r ( -0.7391304347826086 );

    CPPUNIT_ASSERT_EQUAL ( -17ll, r.nominator() );
    CPPUNIT_ASSERT_EQUAL ( 23ll, r.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 77227930ll, m_sqrt2.nominator() );
    CPPUNIT_ASSERT_EQUAL ( 54608393ll, m_sqrt2.denominator() );

    Rational<long long> pi ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( 245850922ll, pi.nominator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779ll, pi.denominator() );

    CPPUNIT_ASSERT_EQUAL ( M_PI, static_cast<double> ( pi ) );
}

void RationalTest::testAssignedFromDouble() {

    Rational<long long> p = 19.0/51.0;

    CPPUNIT_ASSERT_EQUAL ( 19ll, p.nominator() );
    CPPUNIT_ASSERT_EQUAL ( 51ll, p.denominator() );

    Rational<long long> q = 516901.0/740785.0;

    CPPUNIT_ASSERT_EQUAL ( 516901ll, q.nominator() );
    CPPUNIT_ASSERT_EQUAL ( 740785ll, q.denominator() );

    Rational<long long> r = -0.7391304347826086;

    CPPUNIT_ASSERT_EQUAL ( -17ll, r.nominator() );
    CPPUNIT_ASSERT_EQUAL ( 23ll, r.denominator() );

    Rational<long long> pi = M_PI;

    CPPUNIT_ASSERT_EQUAL ( 245850922ll, pi.nominator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779ll, pi.denominator() );
}

void RationalTest::testAddition() {

    const Rational<long long> a ( 17, 21 );
    const Rational<long long> b ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( 31ll, ( a + b ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( 15ll, ( a + b ).denominator() );

    const Rational<long long> c ( 1, 6 );
    const Rational<long long> d ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( 3ll, ( c + d ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( 10ll, ( c + d ).denominator() );
}

void RationalTest::testSubtraction() {

    const Rational<long long> a ( 17, 21 );
    const Rational<long long> b ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( 47ll, ( a - b ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( -105ll, ( a - b ).denominator() );

    const Rational<long long> c ( 1, 6 );
    const Rational<long long> d ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( 1ll, ( c - d ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( 30ll, ( c - d ).denominator() );
}

void RationalTest::testMultiplication() {

    const Rational<long long> a ( 2, 8 );
    const Rational<long long> b ( 7, 3 );

    CPPUNIT_ASSERT_EQUAL ( 7ll, ( a * b ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( 12ll, ( a * b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( 2.0f, static_cast<float> ( m_sqrt2 * m_sqrt2 ) );
}

void RationalTest::testDivision() {

    const Rational<long long> a ( 2, 8 );
    const Rational<long long> b ( 7, 3 );

    CPPUNIT_ASSERT_EQUAL ( 3ll, ( a / b ).nominator() );
    CPPUNIT_ASSERT_EQUAL ( 28ll, ( a / b ).denominator() );
}

void RationalTest::testRelOps() {

    const Rational<long long> a ( 1, 4 );
    const Rational<long long> b ( 1, 2 );

    CPPUNIT_ASSERT ( a < b );
    CPPUNIT_ASSERT ( a <= b );

    CPPUNIT_ASSERT ( b > a );
    CPPUNIT_ASSERT ( b >= a );

    const Rational<long long> c ( 2, 4 );

    CPPUNIT_ASSERT ( c == b );
    CPPUNIT_ASSERT ( b == c );

    CPPUNIT_ASSERT ( a != b );
    CPPUNIT_ASSERT ( b != a );

    CPPUNIT_ASSERT ( b <= c );
    CPPUNIT_ASSERT ( c <= b );
    CPPUNIT_ASSERT ( b >= c );
    CPPUNIT_ASSERT ( c >= b );

    const Rational<long long> d ( 2, 4 );
    const Rational<long long> e ( 2, -4 );

    CPPUNIT_ASSERT ( d > e );
    CPPUNIT_ASSERT ( e < d );

    const Rational<long long> f ( -2, 4 );

    CPPUNIT_ASSERT ( f == e );
    CPPUNIT_ASSERT ( f >= e );
    CPPUNIT_ASSERT ( f <= e );

    CPPUNIT_ASSERT ( e == f );
    CPPUNIT_ASSERT ( e >= f );
    CPPUNIT_ASSERT ( e <= f );

    const Rational<long long> g ( -3, 4 );

    CPPUNIT_ASSERT ( g < d );
    CPPUNIT_ASSERT ( d > g );
}

void RationalTest::testGlobalOps() {

    double a = 0.5;
    a += Rational<long long> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 1.0, a );

    double b = a + Rational<long long> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 1.0, a );
    CPPUNIT_ASSERT_EQUAL ( 1.5, b );

    a -= Rational<long long> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.5, a );

    b = a - Rational<long long> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.5, a );
    CPPUNIT_ASSERT_EQUAL ( 0.0, b );

    a *= Rational<long long> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.25, a );

    b = a * Rational<long long> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.25, a );
    CPPUNIT_ASSERT_EQUAL ( 0.125, b );

    a /= Rational<long long> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.5, a );

    b = a / Rational<long long> ( 1, 2 );

    CPPUNIT_ASSERT_EQUAL ( 0.5, a );
    CPPUNIT_ASSERT_EQUAL ( 1.0, b );

    CPPUNIT_ASSERT ( 0.5 == Rational<long long> ( 1, 2 ) );
    CPPUNIT_ASSERT ( Rational<long long> ( 1, 2 ) == 0.5 );

    CPPUNIT_ASSERT ( 0.5 != Rational<long long> ( 11, 23 ) );
    CPPUNIT_ASSERT ( Rational<long long> ( 11, 23 ) != 0.5 );

    CPPUNIT_ASSERT ( 0.25 < Rational<long long> ( 1, 2 ) );
    CPPUNIT_ASSERT ( Rational<long long> ( 1, 2 ) > 0.25 );

    CPPUNIT_ASSERT ( 0.5 >= Rational<long long> ( 1, 2 ) );
    CPPUNIT_ASSERT ( Rational<long long> ( 1, 2 ) <= 0.5 );

    CPPUNIT_ASSERT ( 0.25 <= Rational<long long> ( 1, 2 ) );
    CPPUNIT_ASSERT ( Rational<long long> ( 1, 2 ) >= 0.25 );
}

void RationalTest::testIOStreamOps() {

    std::ostringstream os;
    os << Rational<long long> ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "245850922/78256779" ), os.str() );

    os.str ( "" );
    os << Rational<unsigned long> ( 280.0f/375.0f );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "56/75" ), os.str() );

    std::istringstream is ( "3.14159265358979323846" );
    Rational<long long> in_pi;

    is >> in_pi;

    CPPUNIT_ASSERT_EQUAL ( 245850922ll, in_pi.nominator() );
    CPPUNIT_ASSERT_EQUAL ( 78256779ll, in_pi.denominator() );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
