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

#include <cppunit/extensions/HelperMacros.h>

#include "rational.h"

class RationalTest : public CppUnit::TestFixture
{
    CPPUNIT_TEST_SUITE ( RationalTest );
    CPPUNIT_TEST ( testNullRational );
    CPPUNIT_TEST ( testConstruct );
    CPPUNIT_TEST ( testConstructFromDouble );
    CPPUNIT_TEST ( testAssignedFromDouble );
    CPPUNIT_TEST ( testAddition );
    CPPUNIT_TEST ( testSubtraction );
    CPPUNIT_TEST ( testMultiplication );
    CPPUNIT_TEST ( testDivision );
    CPPUNIT_TEST ( testRelOps );
    CPPUNIT_TEST ( testGlobalOps );
    CPPUNIT_TEST ( testIOStreamOps );
    CPPUNIT_TEST_SUITE_END();

public:
    void setUp();
    void tearDown();

    void testNullRational();
    void testConstruct();
    void testConstructFromDouble();
    void testAssignedFromDouble();
    void testAddition();
    void testSubtraction();
    void testMultiplication();
    void testDivision();
    void testRelOps();
    void testGlobalOps();
    void testIOStreamOps();

private:
    Commons::Math::Rational<long long> m_nullRational;
    Commons::Math::Rational<long long> m_sqrt2;
};

#endif /* RATIONALTESTCASE_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
