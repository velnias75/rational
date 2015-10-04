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

#include "gmptest.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( GMPTest );

using namespace Commons::Math;

GMPTest::GMPTest() : CppUnit::TestFixture() {}

void GMPTest::setUp() {
    m_sqrt2 = unchecked_sqrt ( mpf_class ( std::sqrt ( 2.0 ) ) );
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
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
