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

#ifndef GMPTESTCASE_H
#define GMPTESTCASE_H

#include <cppunit/extensions/HelperMacros.h>

#include "gmp_rational.h"

#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
class GMPTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( GMPTest );
    CPPUNIT_TEST ( testConstruct );
    CPPUNIT_TEST ( testConstructFromDouble );
    CPPUNIT_TEST ( testConstructFrom_mpf_class );
    CPPUNIT_TEST ( testAddition );
    CPPUNIT_TEST ( testSubtraction );
    CPPUNIT_TEST ( testMultiplication );
    CPPUNIT_TEST ( testInvert );
    CPPUNIT_TEST ( testDivision );
    CPPUNIT_TEST ( testModulo );
    CPPUNIT_TEST ( testIncDec );
    CPPUNIT_TEST ( testRelOps );
    CPPUNIT_TEST ( testString );
    CPPUNIT_TEST ( testIOStreamOps );
    CPPUNIT_TEST ( testAlgorithm );
    CPPUNIT_TEST ( testStdMath );
    CPPUNIT_TEST ( testGoldenRatio );
    CPPUNIT_TEST_SUITE_END();

public:
    typedef mpz_class rational_type;

    GMPTest();

    void setUp();
    void tearDown();

    void testConstruct();
    void testConstructFromDouble();
    void testConstructFrom_mpf_class();
    void testAddition();
    void testSubtraction();
    void testMultiplication();
    void testInvert();
    void testDivision();
    void testModulo();
    void testIncDec();
    void testRelOps();
    void testString();
    void testIOStreamOps();
    void testAlgorithm();
    void testStdMath();
    void testGoldenRatio();

private:
    typedef Commons::Math::Rational<rational_type, Commons::Math::GCD_stein> unchecked_sqrt;
    typedef std::vector<unchecked_sqrt> rat_vector_sqrt;
    typedef std::vector<Commons::Math::gmp_rational> rat_vector;

    unchecked_sqrt m_sqrt2;
    rat_vector_sqrt m_twosqrt;
    rat_vector m_onethird;
    rat_vector m_oneseventh;
};
#pragma GCC diagnostic pop

#endif /* GMPTESTCASE_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
