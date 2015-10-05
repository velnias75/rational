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

    m_sqrt2 = unchecked_sqrt ( mpf_class ( std::sqrt ( 2.0 ) ) );

    std::fill_n ( std::back_inserter ( m_twosqrt ), 2, m_sqrt2 );
    std::fill_n ( std::back_inserter ( m_onethird ), 3, rat_vector::value_type ( 1, 3 ) );
    std::fill_n ( std::back_inserter ( m_oneseventh ), 7, rat_vector::value_type ( 1, 7 ) );
}

void GMPTest::tearDown() {}

void GMPTest::testConstruct() {

#ifdef __EXCEPTIONS
    CPPUNIT_ASSERT_THROW ( Rational<rational_type> r ( 1, 0 ), std::domain_error );
#endif

    CPPUNIT_ASSERT_EQUAL (
        std::string ( "6933950472684408894045574639154274115485852426842518240159" ),
        m_sqrt2.numerator().get_str() );

    CPPUNIT_ASSERT_EQUAL (
        std::string ( "4903043399646811761765323314534746386982093618212104008044" ),
        m_sqrt2.denominator().get_str() );

    const mpf_class a = Rational<rational_type> ( 1, 2 );
    const mpf_class b = Rational<rational_type> ( 1, -2 );
    const mpf_class c = Rational<rational_type> ( -1, 2 );
    const mpf_class d = Rational<rational_type> ( -1, -2 );

    CPPUNIT_ASSERT_EQUAL ( 0.5, a.get_d() );
    CPPUNIT_ASSERT_EQUAL ( -0.5, b.get_d() );
    CPPUNIT_ASSERT_EQUAL ( -0.5, c.get_d() );
    CPPUNIT_ASSERT_EQUAL ( 0.5, d.get_d() );

    CPPUNIT_ASSERT_EQUAL ( -3l, Rational<rational_type> ( 6, -8 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4l, Rational<rational_type> ( 6, -8 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 7l, Rational<rational_type> ( 14, 24 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, Rational<rational_type> ( 14, 24 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 7l, Rational<rational_type> ( 2, 1, 3 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 3l, Rational<rational_type> ( 2, 1, 3 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 86l, Rational<rational_type> ( 18, 4, -5 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 5l, Rational<rational_type> ( 18, 4, -5 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 86l, Rational<rational_type> ( 18, -4, 5 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 5l, Rational<rational_type> ( 18, -4, 5 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( -86l, Rational<rational_type> ( -18, 4, 5 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 5l, Rational<rational_type> ( -18, 4, 5 ).denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( -94l, Rational<rational_type> ( -18, 4, -5 ).numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 5l, Rational<rational_type> ( -18, 4, -5 ).denominator().get_si() );
}

void GMPTest::testConstructFromDouble() {

    const Rational<rational_type, GCD_euclid> p ( 19.0 / 51.0 );

    CPPUNIT_ASSERT_EQUAL ( 19l, p.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 51l, p.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> q ( 516901.0 / 740785.0 );

    CPPUNIT_ASSERT_EQUAL ( 516901l, q.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 740785l, q.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> r ( -0.7391304347826086 );

    CPPUNIT_ASSERT_EQUAL ( -17l, r.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 23l, r.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> s ( 0.0 );

    CPPUNIT_ASSERT_EQUAL ( 0l, s.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, s.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> pi ( M_PI );

    CPPUNIT_ASSERT_EQUAL ( 245850922l, pi.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 78256779l, pi.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> t ( 1.0 );

    CPPUNIT_ASSERT_EQUAL ( 1l, t.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, t.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> u ( 2.0 );

    CPPUNIT_ASSERT_EQUAL ( 2l, u.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, u.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> v ( -8 );

    CPPUNIT_ASSERT_EQUAL ( -8l, v.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, v.denominator().get_si() );
}

void GMPTest::testConstructFrom_mpf_class() {

    const Rational<rational_type, GCD_euclid> p ( mpf_class ( mpf_class ( 19.0 ) /
            mpf_class ( 51.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 19l, p.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 51l, p.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> q ( mpf_class ( mpf_class ( 516901.0 ) /
            mpf_class ( 740785.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 516901l, q.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 740785l, q.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> r ( mpf_class ( -0.7391304347826086, 65 ) );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-665749510133023" ), r.numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "900719925474090" ), r.denominator().get_str() );

    const Rational<rational_type, GCD_euclid> s ( mpf_class ( 0.0 ) );

    CPPUNIT_ASSERT_EQUAL ( 0l, s.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, s.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> pi ( mpf_class ( M_PI ) );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "9978066541" ), pi.numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "3176117225" ), pi.denominator().get_str() );

    const Rational<rational_type, GCD_euclid> t ( mpf_class ( 1.0 ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, t.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, t.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> u ( mpf_class ( 2.0 ) );

    CPPUNIT_ASSERT_EQUAL ( 2l, u.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, u.denominator().get_si() );

    const Rational<rational_type, GCD_euclid> v ( mpf_class ( -8 ) );

    CPPUNIT_ASSERT_EQUAL ( -8l, v.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, v.denominator().get_si() );
}

void GMPTest::testAlgorithm() {

    const mpf_class r =
        std::accumulate ( m_twosqrt.begin(), m_twosqrt.end(),
                          unchecked_sqrt ( 1, 1 ), std::multiplies<unchecked_sqrt > () );
    mp_exp_t exp;

    CPPUNIT_ASSERT_EQUAL ( std::string ( "2" ), r.get_str ( exp, 10, 4 ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           Rational<rat_vector::value_type::integer_type>(),
                           std::plus<Rational<rat_vector::value_type::integer_type> >() )
                           .numerator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           Rational<rat_vector::value_type::integer_type>(),
                           std::plus<Rational<rat_vector::value_type::integer_type> >() )
                           .denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_onethird.begin(), m_onethird.end(),
                           Rational<rat_vector::value_type::integer_type>(),
                           std::plus<Rational<rat_vector::value_type::integer_type> >() )
                           .numerator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 1l, std::accumulate ( m_oneseventh.begin(), m_oneseventh.end(),
                           Rational<rat_vector::value_type::integer_type>(),
                           std::plus<Rational<rat_vector::value_type::integer_type> >() )
                           .denominator().get_si() );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
