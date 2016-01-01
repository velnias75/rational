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

#include <cln/float_io.h>

#include "clntest.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( CLNTest );

using namespace Commons::Math;

CLNTest::CLNTest()
    : CppUnit::TestFixture(), m_sqrt2(), m_twosqrt(), m_onethird(), m_oneseventh() {}

void CLNTest::setUp() {

    m_sqrt2 = unchecked_sqrt ( cln::sqrt ( "2.0L+0_6" ) );

    std::fill_n ( std::back_inserter ( m_twosqrt ), 2, m_sqrt2 );
    std::fill_n ( std::back_inserter ( m_onethird ), 3, rat_vector::value_type ( 1, 3 ) );
    std::fill_n ( std::back_inserter ( m_oneseventh ), 7, rat_vector::value_type ( 1, 7 ) );
}

void CLNTest::tearDown() {}

void CLNTest::testConstruct() {

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( cln_rational r ( 1, 0 ), std::domain_error );
#endif

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "131836323" ), m_sqrt2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "93222358" ), m_sqrt2.denominator() );

    const cln::cl_F &a ( cln_rational ( 1, 2 ) );
    const cln::cl_F &b ( cln_rational ( 1, -2 ) );
    const cln::cl_F &c ( cln_rational ( -1, 2 ) );
    const cln::cl_F &d ( cln_rational ( -1, -2 ) );

    CPPUNIT_ASSERT_EQUAL ( 0.5, cln::double_approx ( a ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, cln::double_approx ( b ) );
    CPPUNIT_ASSERT_EQUAL ( -0.5, cln::double_approx ( c ) );
    CPPUNIT_ASSERT_EQUAL ( 0.5, cln::double_approx ( d ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -3l ), cln_rational ( 6, -8 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 4l ), cln_rational ( 6, -8 ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 7l ), cln_rational ( 14, 24 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 12l ), cln_rational ( 14, 24 ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 7l ), cln_rational ( 2, 1, 3 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3l ), cln_rational ( 2, 1, 3 ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 86l ), cln_rational ( 18, 4, -5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 5l ), cln_rational ( 18, 4, -5 ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 86l ), cln_rational ( 18, -4, 5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 5l ), cln_rational ( 18, -4, 5 ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -86l ), cln_rational ( -18, 4, 5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 5l ), cln_rational ( -18, 4, 5 ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -94l ), cln_rational ( -18, 4, -5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 5l ), cln_rational ( -18, 4, -5 ).denominator() );
}

void CLNTest::testConstructFromDouble() {

    const cln_rational &p ( 19.0 / 51.0 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 19l ), p.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 51l ), p.denominator() );

    const cln_rational &q ( 516901.0 / 740785.0 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 516901l ), q.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 740785l ), q.denominator() );

    const cln_rational &r ( -0.7391304347826086 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -17l ), r.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 23l ), r.denominator() );

    const cln_rational &s ( 0.0 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 0l ), s.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), s.denominator() );

    const cln_rational &pi ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 245850922l ), pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 78256779l ), pi.denominator() );

    const cln_rational &t ( 1.0 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), t.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), t.denominator() );

    const cln_rational &u ( 2.0 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), u.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), u.denominator() );

    const cln_rational &v ( -8 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -8l ), v.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), v.denominator() );
}

void CLNTest::testConstructFrom_cl_F_class() {

    const cln_rational &o ( cln::cl_F ( "0.33333333333333333L0_16" ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), o.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3l ), o.denominator() );

    const cln_rational &p ( cln::cl_F ( 19.0 ) / cln::cl_F ( 51.0 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 19l ), p.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 51l ), p.denominator() );

    const cln_rational &q ( cln::cl_F ( 516901.0 ) / cln::cl_F ( 740785.0 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 516901l ), q.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 740785l ), q.denominator() );

    const cln_rational &r ( cln::cl_F ( "-0.7391304347826086L+0_65" ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -17 ), r.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 23 ), r.denominator() );

    const cln_rational &s ( cln::cl_F ( 0.0 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 0l ), s.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), s.denominator() );

    const cln_rational &pi ( cln::cl_F ( M_PI ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "245850922" ), pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "78256779" ), pi.denominator() );

    const cln_rational &t ( cln::cl_F ( 1.0 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), t.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), t.denominator() );

    const cln_rational &u ( cln::cl_F ( 2.0 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), u.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), u.denominator() );

    const cln_rational &v ( cln::cl_I ( -8 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -8l ), v.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), v.denominator() );

}

void CLNTest::testAddition() {

    const cln_rational a ( 17, 21 );
    const cln_rational b ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 31l ), ( a + b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 15l ), ( a + b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 31l ), ( b + a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 15l ), ( b + a ).denominator() );

    const cln_rational c ( 1, 6 );
    const cln_rational d ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3l ), ( c + d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 10l ), ( c + d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3l ), ( d + c ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 10l ), ( d + c ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), ( +d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 15l ), ( +d ).denominator() );

    const cln_rational knuth_a ( 7, 66 );
    const cln_rational knuth_b ( 17, 12 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 67l ), ( knuth_a + knuth_b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 44l ), ( knuth_a + knuth_b ).denominator() );
}

void CLNTest::testSubtraction() {

    const cln_rational a ( 17, 21 );
    const cln_rational b ( 44, 35 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -47l ), ( a - b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 105l ), ( a - b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 0l ), ( a - a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), ( a - a ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 47l ), ( b - a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 105l ), ( b - a ).denominator() );

    const cln_rational c ( 1, 6 );
    const cln_rational d ( 2, 15 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), ( c - d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 30l ), ( c - d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -1l ), ( d - c ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 30l ), ( d - c ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -2l ), ( -d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 15l ), ( -d ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), ( d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 15l ), ( d ).denominator() );

}

void CLNTest::testMultiplication() {

    const cln_rational a ( 2, 8 );
    const cln_rational b ( 7, 3 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 7l ), ( a * b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 12l ), ( a * b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 7l ), ( b * a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 12l ), ( b * a ).denominator() );

}

void CLNTest::testInvert() {
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 7l ), cln_rational ( 161, 49 ).invert().numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 23l ), cln_rational ( 161, 49 ).invert().denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -7l ), cln_rational ( 161, -49 ).inverse().numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 23l ), cln_rational ( -161, 49 ).inverse().denominator() );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( cln_rational ().invert(), std::domain_error );
    CPPUNIT_ASSERT_THROW ( cln_rational ().inverse(), std::domain_error );
#endif
}

void CLNTest::testDivision() {

    const cln_rational a ( 2, 8 );
    const cln_rational b ( 7, 3 );
    const cln_rational c ( 0, 1 );
    const cln_rational d ( -7, -3 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3l ), ( a / b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 28l ), ( a / b ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 28l ), ( b / a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3l ), ( b / a ).denominator() );

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( a / c, std::domain_error );
    CPPUNIT_ASSERT_THROW ( a / ( b - d ), std::domain_error );
#endif
}

void CLNTest::testModulo() {

    cln_rational a ( 8, 1 );

    a %= cln_rational ( 3, 1 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), ( a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), ( a ).denominator() );

    cln_rational c ( 41, 7 );
    c %= cln_rational ( 3, 2 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 19l ), ( c ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 14l ), ( c ).denominator() );

    cln_rational d ( 542, 84 );
    cln_rational e ( -65, 28 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -43l ), ( d % e ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 84l ), ( d % e ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 347l ), ( e % d ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 84l ), ( e % d ).denominator() );

    cln_rational h ( 11, 4 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), h.mod().first );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3l ), h.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 4l ), h.mod().second.denominator() );

    cln_rational i ( 11, -4 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -2l ), i.mod().first );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -3l ), i.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 4l ), i.mod().second.denominator() );

    cln_rational j ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), j.mod().first );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), j.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 4l ), j.mod().second.denominator() );

    cln_rational k ( -18, 8 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -2l ), k.mod().first );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -1l ), k.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 4l ), k.mod().second.denominator() );

    cln_rational l ( 1, 8 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 0l ), l.mod().first );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), l.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 8l ), l.mod().second.denominator() );
}

void CLNTest::testIncDec() {

    cln_rational a ( 2, 4 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3l ), ( ++a ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), ( a++ ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 5l ), a.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), a.denominator() );

    cln_rational b ( 2, 4 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -1l ), ( --b ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), ( b-- ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -3l ), b.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), b.denominator() );
}

void CLNTest::testRelOps() {

    const cln_rational a ( 1, 4 );
    const cln_rational b ( 1, 2 );

    CPPUNIT_ASSERT ( a < b );
    CPPUNIT_ASSERT ( a <= b );

    CPPUNIT_ASSERT ( b > a );
    CPPUNIT_ASSERT ( b >= a );

    const cln_rational c ( 2, 4 );

    CPPUNIT_ASSERT ( c == b );
    CPPUNIT_ASSERT ( b == c );

    CPPUNIT_ASSERT ( a != b );
    CPPUNIT_ASSERT ( b != a );

    CPPUNIT_ASSERT ( b <= c );
    CPPUNIT_ASSERT ( c <= b );
    CPPUNIT_ASSERT ( b >= c );
    CPPUNIT_ASSERT ( c >= b );

    const cln_rational d ( 2, 4 );
    const cln_rational e ( 2, -4 );

    CPPUNIT_ASSERT ( d > e );
    CPPUNIT_ASSERT ( e < d );

    const cln_rational f ( -2, 4 );

    CPPUNIT_ASSERT ( f == e );
    CPPUNIT_ASSERT ( f >= e );
    CPPUNIT_ASSERT ( f <= e );

    CPPUNIT_ASSERT ( e == f );
    CPPUNIT_ASSERT ( e >= f );
    CPPUNIT_ASSERT ( e <= f );

    const cln_rational g ( -3, 4 );

    CPPUNIT_ASSERT ( g < d );
    CPPUNIT_ASSERT ( d > g );

}

void CLNTest::testString() {

    std::ostringstream os_f, os_m;

    const cln_rational h ( 11, 4 );

    os_f << h.str();
    os_m << h.str ( true );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/4" ), os_f.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "2 3/4" ), os_m.str() );

    const cln_rational i ( 11, -4 );

    os_f.str ( "" );
    os_m.str ( "" );

    os_f << i.str();
    os_m << i.str ( true );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-11/4" ), os_f.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2 3/4" ), os_m.str ( ) );

    const cln_rational j ( 18, 8 );

    os_f.str ( "" );
    os_m.str ( "" );

    os_f << j.str();
    os_m << j.str ( true );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "9/4" ), os_f.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "2 1/4" ), os_m.str ( ) );

    const cln_rational k ( -18, 8 );

    os_f.str ( "" );
    os_m.str ( "" );

    os_f << k.str();
    os_m << k.str ( true );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-9/4" ), os_f.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2 1/4" ), os_m.str ( ) );

    const cln_rational l ( 1, 8 );

    os_f.str ( "" );
    os_m.str ( "" );

    os_f << l.str();
    os_m << l.str ( true );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "1/8" ), os_f.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "1/8" ), os_m.str() );

    const cln_rational m ( 8, 1 );

    os_f.str ( "" );
    os_m.str ( "" );

    os_f << m.str();
    os_m << m.str ( true );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "8" ), os_f.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "8" ), os_m.str() );

    const cln_rational n ( 8, 2, 1 );

    os_f.str ( "" );
    os_m.str ( "" );

    os_f << n.str();
    os_m << n.str ( true );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "10" ), os_f.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "10" ), os_m.str() );

}

void CLNTest::testIOStreamOps() {

    std::istringstream real_in ( "0.33333333333333333" );
    cln_rational real_rat;

    real_in >> real_rat;

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "1" ), real_rat.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "3" ), real_rat.denominator() );

    std::ostringstream os;
    os << cln_rational ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "245850922/78256779" ), os.str() );

    os.str ( "" );
    os << cln_rational ( 280.0f/375.0f );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "56/75" ), os.str() );

    cln_rational in_pi;
    std::istringstream is ( "6.14159265358979323846 - (1 + 2)" );

    is >> in_pi;

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "245850922" ), in_pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "78256779" ), in_pi.denominator() );

}

#pragma GCC diagnostic ignored "-Wuseless-cast"
#pragma GCC diagnostic push
void CLNTest::testAlgorithm() {

    std::ostringstream os;

    const cln::cl_F &r (
        std::accumulate ( m_twosqrt.begin(), m_twosqrt.end(),
                          unchecked_sqrt ( 1, 1 ), std::multiplies<unchecked_sqrt> () ) );

    os << r;

    CPPUNIT_ASSERT_EQUAL ( std::string ( "2.00000000000000011506939563983927016233L0" ), os.str() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           rat_vector::value_type(), std::plus<rat_vector::value_type>() )
                           .numerator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           rat_vector::value_type(), std::plus<rat_vector::value_type>() )
                           .denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           rat_vector::value_type(), std::plus<rat_vector::value_type>() )
                           .numerator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), std::accumulate ( m_oneseventh.begin(),
                           m_oneseventh.end(), rat_vector::value_type(),
                           std::plus<rat_vector::value_type>() ).denominator() );

    const cln::cl_I cf_pi[] = {  3,  7, 15,  1, 292,  1, 1,   1,  2,  1,  3, 1, 14, 2, 1,  1, 2, 2,
                                 2,  2,  1, 84,   2,  1, 1,  15,  3, 13,  1, 4,  2, 6, 6, 99, 1, 2,
                                 2,  6,  3,  5,   1,  1, 6,   8,  1,  7,  1, 2,  3, 7, 1,  2, 1, 1,
                                 12, 1,  1,  1,   3,  1, 1,   8,  1,  1,  2, 1,  6, 1, 1,  5, 2, 2,
                                 3,  1,  2,  4,   4, 16, 1, 161, 45,  1, 22, 1,  2, 2, 1,  4, 1, 2,
                                 24, 1,  2,  1,   3,  1, 3
                              };

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "227159758552934520439668309319746303422708645581861" ),
                           cf ( cf_pi, cf_pi + 97 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "72307196890521956737416455481060519150048966236850" ),
                           cf ( cf_pi, cf_pi + 97 ).denominator() );

    std::vector<cln_rational::integer_type> o_pi;
    seq ( cf ( cf_pi, cf_pi + 97 ), std::back_inserter ( o_pi ) );

    CPPUNIT_ASSERT_EQUAL ( static_cast<std::vector<cln_rational::integer_type>::size_type> ( 97u ), o_pi.size() );
    CPPUNIT_ASSERT ( std::equal ( o_pi.begin(), o_pi.end(), cf_pi ) );
}
#pragma GCC diagnostic pop

void CLNTest::testStdMath() {

    rational_type rt;

    CPPUNIT_ASSERT_EQUAL ( std::string ( "2/3" ), std::modf ( cln_rational ( 11, 3 ) ,
                           &rt ).str() );
    CPPUNIT_ASSERT_EQUAL ( rational_type ( 3 ), rt );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), cln_rational ( 11, -3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), cln_rational ( -11, 3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), cln_rational ( 11, 3 ).abs().str() );

    const cln_rational &a ( cln_rational::rf_info ( 142857 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1 ), a.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 7 ), a.denominator() );

    const cln_rational &b ( cln_rational::rf_info ( 34 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 34 ), b.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 99 ), b.denominator() );

    const cln_rational &c ( cln_rational::rf_info ( 123456789 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 13717421 ), c.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 111111111 ), c.denominator() );

    const cln_rational &d ( cln_rational::rf_info ( 12, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 4 ), d.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 333 ), d.denominator() );

    const cln_rational &ex ( cln_rational::rf_info ( 6, 0, 1111 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 667 ), ex.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 6000 ), ex.denominator() );

    const cln_rational &f ( cln_rational::rf_info ( 1, 2, 3, 4 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1499 ), f.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 49950000 ), f.denominator() );

    const cln_rational &g ( cln_rational::rf_info ( 6, 0, 0, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1 ), g.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 15 ), g.denominator() );

    const cln_rational &h ( cln_rational::rf_info ( 6, 0, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1 ), h.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 6 ), h.denominator() );

    const cln_rational &i ( cln_rational::rf_info ( 1, 1 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1 ), i.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 99 ), i.denominator() );

    const cln_rational &j ( cln_rational::rf_info ( 1 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1 ), j.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 9 ), j.denominator() );

    cln_rational::rf_info dc;

    const cln_rational k ( 7, 13 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 0 ), k.decompose ( dc ) );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 7 ), cln_rational ( dc ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 13 ), cln_rational ( dc ).denominator() );

    const cln_rational l ( 88, 100 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 0 ), l.decompose ( dc ) );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 22 ), cln_rational ( dc ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 25 ), cln_rational ( dc ).denominator() );

    const cln_rational m ( 8, 3 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2 ), m.decompose ( dc ) );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2 ), cln_rational ( dc ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3 ), cln_rational ( dc ).denominator() );

    const cln_rational n ( "(70/2) - (1741832/249975)" );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 28 ), n.decompose ( dc ) );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3 ), dc.pre );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 1 ), dc.pre_leading_zeros );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1975 ), dc.reptend );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 0 ), dc.leading_zeros );

    const cln_rational q ( "123.32 / (12453/370)" );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 228142 ), q.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 62265 ), q.denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3 ), q.decompose ( dc ) );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 6 ), dc.pre );

    std::ostringstream os;

    os << dc.reptend;

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
                                         "515217216734923311651810808" ), os.str() );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 0 ), dc.pre_leading_zeros );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 1776 ), dc.reptend_digits.size() );
    CPPUNIT_ASSERT_EQUAL ( std::size_t ( 0 ), dc.leading_zeros );

    const cln_rational s ( 3, 4 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 81l ), s.pow ( 4 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 256l ), s.pow ( 4 ).denominator() );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 243l ), s.pow ( 5 ).numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1024l ), s.pow ( 5 ).denominator() );

#ifdef __EXCEPTIONS
    const cln_rational t ( 3, 4 );

    CPPUNIT_ASSERT_THROW ( t.pow ( 0 ), std::domain_error );
    CPPUNIT_ASSERT_THROW ( t.pow ( -8 ), std::domain_error );
#endif

    const cln_rational u ( 2, 1 );

    os.str ( "" );
    os << u.sqrt().numerator();

    CPPUNIT_ASSERT_EQUAL ( std::string ( "4946041176255201878775086487573351061418968498177" ),
                           os.str() );

    os.str ( "" );
    os << u.sqrt().denominator();

    CPPUNIT_ASSERT_EQUAL ( std::string ( "3497379255757941172020851852070562919437964212608" ),
                           os.str() );

    const cln_rational v ( 10, 17 );

    os.str ( "" );
    os << v.sqrt().numerator();

    CPPUNIT_ASSERT_EQUAL ( std::string ( "1983567417147843927170789761" ), os.str() );

    os.str ( "" );
    os << v.sqrt().denominator();

    CPPUNIT_ASSERT_EQUAL ( std::string ( "2586255495350365951590026592" ),  os.str() );
}

void CLNTest::testGoldenRatio() {

    Rational<cln_rational::integer_type, GCD_null> phi ( cln_rational::one_, cln_rational::one_ );

    for ( std::size_t i = 0u; i < 1024u; ++i ) ( ++phi ).invert();

    std::ostringstream os;

    os << phi.inverse().numerator();

    CPPUNIT_ASSERT_EQUAL ( std::string ( "1179869281805523255014757888412586560808902854456091" \
                                         "3468519228968187430794620907976123201977895385245239" \
                                         "7050828306569046301783141598663704952115390234610526" \
                                         "8281123032179655593090772272438413164852733945840731" \
                                         "7543768" ), os.str() );
    os.str ( "" );
    os << phi.inverse().denominator();

    CPPUNIT_ASSERT_EQUAL ( std::string ( "7291993184377412737043195648396979558721167948342308" \
                                         "6377162058185874001489121865798744093687543548489948" \
                                         "3181625031189341064810479244078947534047137736685242" \
                                         "0526027975140687031196633477605718294523235826853392" \
                                         "138525" ), os.str() );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
