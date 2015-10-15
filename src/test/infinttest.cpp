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

#include "infinttest.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( InfIntTest );

using namespace Commons::Math;

InfIntTest::InfIntTest() : CppUnit::TestFixture(), m_sqrt2(), m_onethird(), m_oneseventh() {}

void InfIntTest::setUp() {

    using namespace std;

    m_sqrt2 = unchecked_sqrt ( std::sqrt ( 2.0l ) );

    std::fill_n ( std::back_inserter ( m_twosqrt ), 2, m_sqrt2 );
    std::fill_n ( std::back_inserter ( m_onethird ), 3, rat_vector::value_type ( 1, 3 ) );
    std::fill_n ( std::back_inserter ( m_oneseventh ), 7, rat_vector::value_type ( 1, 7 ) );
}

void InfIntTest::tearDown() {}

void InfIntTest::testConstruct() {

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( infint_rational r ( 1, 0 ), std::domain_error );
#endif

    CPPUNIT_ASSERT_EQUAL ( InfInt ( 6333631924 ), m_sqrt2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( InfInt ( 4478554083 ), m_sqrt2.denominator() );

    const long double &a ( infint_rational ( 1, 2 ) );
    const long double &b ( infint_rational ( 1, -2 ) );
    const long double &c ( infint_rational ( -1, 2 ) );
    const long double &d ( infint_rational ( -1, -2 ) );

    CPPUNIT_ASSERT_EQUAL ( 0.5l, a );
    CPPUNIT_ASSERT_EQUAL ( -0.5l, b );
    CPPUNIT_ASSERT_EQUAL ( -0.5l, c );
    CPPUNIT_ASSERT_EQUAL ( 0.5l, d );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-3" ), infint_rational ( 6, -8 ).numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 4l, infint_rational ( 6, -8 ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 7l, infint_rational ( 14, 24 ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 12l, infint_rational ( 14, 24 ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 7l, infint_rational ( 2, 1, 3 ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 3l, infint_rational ( 2, 1, 3 ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 86l, infint_rational ( 18, 4, -5 ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 5l, infint_rational ( 18, 4, -5 ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 86l, infint_rational ( 18, -4, 5 ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 5l, infint_rational ( 18, -4, 5 ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-86" ),
                           infint_rational ( -18, 4, 5 ).numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 5l, infint_rational ( -18, 4, 5 ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-94" ),
                           infint_rational ( -18, 4, -5 ).numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 5l, infint_rational ( -18, 4, -5 ).denominator().toLong() );
}

void InfIntTest::testConstructFromDouble() {

    const Rational<rational_type, GCD_euclid> &p ( 19.0 / 51.0 );

    CPPUNIT_ASSERT_EQUAL ( 19l, p.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 51l, p.denominator().toLong() );

    const Rational<rational_type, GCD_euclid> &q ( 516901.0 / 740785.0 );

    CPPUNIT_ASSERT_EQUAL ( 516901l, q.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 740785l, q.denominator().toLong() );

    const Rational<rational_type, GCD_euclid> &r ( -0.7391304347826086 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-17" ), r.numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 23l, r.denominator().toLong() );

    const Rational<rational_type, GCD_euclid> &s ( 0.0 );

    CPPUNIT_ASSERT_EQUAL ( 0l, s.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 1l, s.denominator().toLong() );

    const Rational<rational_type, GCD_euclid> &pi ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( 245850922l, pi.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 78256779l, pi.denominator().toLong() );

    const Rational<rational_type, GCD_euclid> &t ( 1.0 );

    CPPUNIT_ASSERT_EQUAL ( 1l, t.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 1l, t.denominator().toLong() );

    const Rational<rational_type, GCD_euclid> &u ( 2.0 );

    CPPUNIT_ASSERT_EQUAL ( 2l, u.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 1l, u.denominator().toLong() );

    const Rational<rational_type, GCD_euclid> &v ( -8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-8" ), v.numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 1l, v.denominator().toLong() );
}

void InfIntTest::testAddition() {

    const Rational<rational_type, GCD_euclid> a ( 17, 21 );
    const Rational<rational_type, GCD_euclid> b ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( 31l, ( a + b ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 15l, ( a + b ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 31l, ( b + a ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 15l, ( b + a ).denominator().toLong() );

    const Rational<rational_type, GCD_euclid> c ( 1, 6 );
    const Rational<rational_type, GCD_euclid> d ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( 3l, ( c + d ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 10l, ( c + d ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 3l, ( d + c ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 10l, ( d + c ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 2l, ( +d ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 15l, ( +d ).denominator().toLong() );

    const Rational<rational_type, GCD_euclid> knuth_a ( 7, 66 );
    const Rational<rational_type, GCD_euclid> knuth_b ( 17, 12 );

    CPPUNIT_ASSERT_EQUAL ( 67l, ( knuth_a + knuth_b ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 44l, ( knuth_a + knuth_b ).denominator().toLong() );
}

void InfIntTest::testSubtraction() {

    const Rational<rational_type, GCD_euclid> a ( 17, 21 );
    const Rational<rational_type, GCD_euclid> b ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-47" ), ( a - b ).numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 105l, ( a - b ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 0l, ( a - a ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 1l, ( a - a ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 47l, ( b - a ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 105l, ( b - a ).denominator().toLong() );

    const Rational<rational_type, GCD_euclid> c ( 1, 6 );
    const Rational<rational_type, GCD_euclid> d ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( 1l, ( c - d ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 30l, ( c - d ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-1" ), ( d - c ).numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 30l, ( d - c ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2" ), ( -d ).numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 15l, ( -d ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 2l, ( d ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 15l, ( d ).denominator().toLong() );
}

void InfIntTest::testMultiplication() {

    const infint_rational a ( 2, 8 );
    const infint_rational b ( 7, 3 );

    CPPUNIT_ASSERT_EQUAL ( 7l, ( a * b ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 12l, ( a * b ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 7l, ( b * a ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 12l, ( b * a ).denominator().toLong() );
}

void InfIntTest::testInvert() {

    CPPUNIT_ASSERT_EQUAL ( 7l, infint_rational ( 161, 49 ).invert().numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 23l, infint_rational ( 161, 49 ).invert().denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 7l, infint_rational ( 161, 49 ).inverse().numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 23l, infint_rational ( 161, 49 ).inverse().denominator().toLong() );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( infint_rational ().invert(), std::domain_error );
    CPPUNIT_ASSERT_THROW ( infint_rational ().inverse(), std::domain_error );
#endif
}

void InfIntTest::testDivision() {

    const infint_rational a ( 2, 8 );
    const infint_rational b ( 7, 3 );
    const infint_rational c ( 0, 1 );
    const infint_rational d ( -7, -3 );

    CPPUNIT_ASSERT_EQUAL ( 3l, ( a / b ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 28l, ( a / b ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 28l, ( b / a ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 3l, ( b / a ).denominator().toLong() );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( a / c, std::domain_error );
    CPPUNIT_ASSERT_THROW ( a / ( b - d ), std::domain_error );
#endif
}

void InfIntTest::testModulo() {

    infint_rational a ( 8, 1 );

    a %= infint_rational ( 3, 1 );

    CPPUNIT_ASSERT_EQUAL ( 2l, ( a ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 1l, ( a ).denominator().toLong() );

    infint_rational c ( 41, 7 );
    c %= infint_rational ( 3, 2 );

    CPPUNIT_ASSERT_EQUAL ( 19l, ( c ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 14l, ( c ).denominator().toLong() );

    infint_rational d ( 542, 84 );
    infint_rational e ( -65, 28 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-43" ), ( d % e ).numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 84l, ( d % e ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 347l, ( e % d ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 84l, ( e % d ).denominator().toLong() );

    infint_rational h ( 11, 4 );

    CPPUNIT_ASSERT_EQUAL ( 2l, h.mod().first.toLong() );
    CPPUNIT_ASSERT_EQUAL ( 3l, h.mod().second.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 4l, h.mod().second.denominator().toLong() );

    infint_rational i ( 11, -4 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2" ), i.mod().first.toString() );
    CPPUNIT_ASSERT_EQUAL ( 3l, i.mod().second.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 4l, i.mod().second.denominator().toLong() );

    infint_rational j ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( 2l, j.mod().first.toLong() );
    CPPUNIT_ASSERT_EQUAL ( 1l, j.mod().second.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 4l, j.mod().second.denominator().toLong() );

    infint_rational k ( -18, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2" ), k.mod().first.toString() );
    CPPUNIT_ASSERT_EQUAL ( 1l, k.mod().second.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 4l, k.mod().second.denominator().toLong() );

    infint_rational l ( 1, 8 );

    CPPUNIT_ASSERT_EQUAL ( 0l, l.mod().first.toLong() );
    CPPUNIT_ASSERT_EQUAL ( 1l, l.mod().second.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 8l, l.mod().second.denominator().toLong() );
}

void InfIntTest::testIncDec() {

    infint_rational a ( 2, 4 );

    CPPUNIT_ASSERT_EQUAL ( 3l, ( ++a ).numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 2l, ( a++ ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 5l, a.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 2l, a.denominator().toLong() );

    infint_rational b ( 2, 4 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-1" ), ( --b ).numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 2l, ( b-- ).denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-3" ), b.numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 2l, b.denominator().toLong() );
}

void InfIntTest::testRelOps() {

    const infint_rational a ( 1, 4 );
    const infint_rational b ( 1, 2 );

    CPPUNIT_ASSERT ( a < b );
    CPPUNIT_ASSERT ( a <= b );

    CPPUNIT_ASSERT ( b > a );
    CPPUNIT_ASSERT ( b >= a );

    const infint_rational c ( 2, 4 );

    CPPUNIT_ASSERT ( c == b );
    CPPUNIT_ASSERT ( b == c );

    CPPUNIT_ASSERT ( a != b );
    CPPUNIT_ASSERT ( b != a );

    CPPUNIT_ASSERT ( b <= c );
    CPPUNIT_ASSERT ( c <= b );
    CPPUNIT_ASSERT ( b >= c );
    CPPUNIT_ASSERT ( c >= b );

    const infint_rational d ( 2, 4 );
    const infint_rational e ( 2, -4 );

    CPPUNIT_ASSERT ( d > e );
    CPPUNIT_ASSERT ( e < d );

    const infint_rational f ( -2, 4 );

    CPPUNIT_ASSERT ( f == e );
    CPPUNIT_ASSERT ( f >= e );
    CPPUNIT_ASSERT ( f <= e );

    CPPUNIT_ASSERT ( e == f );
    CPPUNIT_ASSERT ( e >= f );
    CPPUNIT_ASSERT ( e <= f );

    const infint_rational g ( -3, 4 );

    CPPUNIT_ASSERT ( g < d );
    CPPUNIT_ASSERT ( d > g );
}

void InfIntTest::testString() {

    const infint_rational h ( 11, 4 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/4" ), h.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "2 3/4" ), h.str ( true ) );

    const infint_rational i ( 11, -4 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-11/4" ), i.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2 3/4" ), i.str ( true ) );

    const infint_rational j ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "9/4" ), j.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "2 1/4" ), j.str ( true ) );

    const infint_rational k ( -18, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-9/4" ), k.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2 1/4" ), k.str ( true ) );

    const infint_rational l ( 1, 8 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "1/8" ), l.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "1/8" ), l.str ( true ) );

    const infint_rational m ( 8, 1 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "8" ), m.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "8" ), m.str ( true ) );

    const infint_rational n ( 8, 2, 1 );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "10" ), n.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "10" ), n.str ( true ) );
}

void InfIntTest::testIOStreamOps() {

    std::ostringstream os;
    os << infint_rational ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "245850922/78256779" ), os.str() );

    os.str ( "" );
    os << infint_rational ( 280.0f/375.0f );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "56/75" ), os.str() );

    std::istringstream is ( "3.14159265358979323846" );
    infint_rational in_pi;

    is >> in_pi;

    CPPUNIT_ASSERT_EQUAL ( 245850922l, in_pi.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 78256779l, in_pi.denominator().toLong() );
}

void InfIntTest::testAlgorithm() {

    const long double &r (
        std::accumulate ( m_twosqrt.begin(), m_twosqrt.end(),
                          unchecked_sqrt ( 1, 1 ), std::multiplies<unchecked_sqrt> () ) );

    CPPUNIT_ASSERT_EQUAL ( 2.0l, r );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           rat_vector::value_type(), std::plus<rat_vector::value_type>() )
                           .numerator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           rat_vector::value_type(), std::plus<rat_vector::value_type>() )
                           .denominator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           rat_vector::value_type(), std::plus<rat_vector::value_type>() )
                           .numerator().toLong() );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_oneseventh.begin(), m_oneseventh.end(),
                           rat_vector::value_type(),std::plus<rat_vector::value_type>() )
                           .denominator().toLong() );
}

void InfIntTest::testStdMath() {

    rational_type rt;

    CPPUNIT_ASSERT_EQUAL ( std::string ( "2/3" ), std::modf ( infint_rational ( 11, 3 ) ,
                           &rt ).str() );
    CPPUNIT_ASSERT_EQUAL ( rational_type ( 3 ), rt );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), infint_rational ( 11, -3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), infint_rational ( -11, 3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), infint_rational ( 11, 3 ).abs().str() );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
