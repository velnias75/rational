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

    // TODO: euclid failed here (modulus-op)

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
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3l ), i.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 4l ), i.mod().second.denominator() );

    cln_rational j ( 18, 8 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2l ), j.mod().first );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), j.mod().second.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 4l ), j.mod().second.denominator() );

    cln_rational k ( -18, 8 );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -2l ), k.mod().first );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), k.mod().second.numerator() );
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
    CPPUNIT_ASSERT_EQUAL ( std::string ( "1/8" ), os_m.str ( ) );

    const cln_rational m ( 8, 1 );

    os_f.str ( "" );
    os_m.str ( "" );

    os_f << m.str();
    os_m << m.str ( true );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "8" ), os_f.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "8" ), os_m.str ( ) );

    const cln_rational n ( 8, 2, 1 );

    os_f.str ( "" );
    os_m.str ( "" );

    os_f << n.str();
    os_m << n.str ( true );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "10" ), os_f.str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "10" ), os_m.str ( ) );

}

void CLNTest::testIOStreamOps() {

    std::ostringstream os;
    os << cln_rational ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "245850922/78256779" ), os.str() );

    os.str ( "" );
    os << cln_rational ( 280.0f/375.0f );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "56/75" ), os.str() );

    cln_rational in_pi;

    ( std::istringstream ( "6.14159265358979323846 - (1 + 2)" ) ) >> in_pi;

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "6167950454" ), in_pi.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "1963319607" ), in_pi.denominator() );

}

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
}

void CLNTest::testStdMath() {

    rational_type rt;

    CPPUNIT_ASSERT_EQUAL ( std::string ( "2/3" ), std::modf ( cln_rational ( 11, 3 ) ,
                           &rt ).str() );
    CPPUNIT_ASSERT_EQUAL ( rational_type ( 3 ), rt );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), cln_rational ( 11, -3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), cln_rational ( -11, 3 ).abs().str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "11/3" ), cln_rational ( 11, 3 ).abs().str() );
}

void CLNTest::testGoldenRatio() {

    const cln_rational one ( cln_rational::one_, cln_rational::one_ );
    cln_rational phi ( one );

    for ( std::size_t i = 0u; i < 1024u; ++i ) ( phi += one ).invert();

    std::ostringstream os;

    os << phi.numerator();

    CPPUNIT_ASSERT_EQUAL ( std::string ( "7291993184377412737043195648396979558721167948342308" \
                                         "6377162058185874001489121865798744093687543548489948" \
                                         "3181625031189341064810479244078947534047137736685242" \
                                         "0526027975140687031196633477605718294523235826853392" \
                                         "138525" ), os.str() );
    os.str ( "" );
    os << phi.denominator();

    CPPUNIT_ASSERT_EQUAL ( std::string ( "1179869281805523255014757888412586560808902854456091" \
                                         "3468519228968187430794620907976123201977895385245239" \
                                         "7050828306569046301783141598663704952115390234610526" \
                                         "8281123032179655593090772272438413164852733945840731" \
                                         "7543768" ), os.str() );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
