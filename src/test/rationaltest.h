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

#ifndef RATIONALTESTCASE_H
#define RATIONALTESTCASE_H

#include <vector>

#include <stdint.h>

#include <cppunit/extensions/HelperMacros.h>

#include "rational.h"

#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
class RationalTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( RationalTest );
    CPPUNIT_TEST ( testNullRational );
    CPPUNIT_TEST ( testConstruct );
    CPPUNIT_TEST ( testConstructFromDouble );
    CPPUNIT_TEST ( testAssignedFromDouble );
    CPPUNIT_TEST ( testAddition );
    CPPUNIT_TEST ( testSubtraction );
    CPPUNIT_TEST ( testMultiplication );
    CPPUNIT_TEST ( testInvert );
    CPPUNIT_TEST ( testDivision );
    CPPUNIT_TEST ( testModulo );
    CPPUNIT_TEST ( testIncDec );
    CPPUNIT_TEST ( testRelOps );
    CPPUNIT_TEST ( testGlobalOps );
    CPPUNIT_TEST ( testIOStreamOps );
    CPPUNIT_TEST ( testPrecision );
    CPPUNIT_TEST ( testAlgorithm );
    CPPUNIT_TEST_SUITE_END();

public:
    typedef int32_t rational_type;

    RationalTest();

    void setUp();
    void tearDown();

    void testNullRational();
    void testConstruct();
    void testConstructFromDouble();
    void testAssignedFromDouble();
    void testAddition();
    void testSubtraction();
    void testMultiplication();
    void testInvert();
    void testDivision();
    void testModulo();
    void testIncDec();
    void testRelOps();
    void testGlobalOps();
    void testIOStreamOps();
    void testPrecision();
    void testAlgorithm();

private:
    Commons::Math::Rational<rational_type> m_nullRational;
    Commons::Math::Rational<uint64_t> m_sqrt2;

    typedef std::vector<Commons::Math::Rational<rational_type> > rat_vector;
    rat_vector m_accu;
    rat_vector m_onethird;

    typedef std::vector<Commons::Math::Rational<uint32_t> > rat_vector_ul;
    rat_vector_ul m_accu_ul;
};
#pragma GCC diagnostic pop

#endif /* RATIONALTESTCASE_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
