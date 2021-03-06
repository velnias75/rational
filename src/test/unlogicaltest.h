/*
 * Copyright 2018 by Heiko Schäfer <heiko@rangun.de>
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

#ifndef UNLOGICALTESTCASE_H
#define UNLOGICALTESTCASE_H

#include <cppunit/extensions/HelperMacros.h>

#include "unlogical_rational.h"

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
#define RATIONAL_OVERRIDE override
#define RATIONAL_FINAL final
#else
#define RATIONAL_OVERRIDE
#define RATIONAL_FINAL
#endif

#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
class UnlogicalTest RATIONAL_FINAL : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( UnlogicalTest );
    CPPUNIT_TEST ( testConstruct );
    CPPUNIT_TEST ( testConstructFromDouble );
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
    typedef Commons::Lab::unlogical::Number<int64_t, 64> rational_type;

    UnlogicalTest();

    void setUp() RATIONAL_OVERRIDE;
    void tearDown() RATIONAL_OVERRIDE;

    void testConstruct();
    void testConstructFromDouble();
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
    typedef Commons::Math::unlogical_rational unchecked_sqrt;
    typedef std::vector<unchecked_sqrt> rat_vector_sqrt;
    typedef std::vector<Commons::Math::unlogical_rational> rat_vector;

    unchecked_sqrt m_sqrt2;
    rat_vector_sqrt m_twosqrt;
    rat_vector m_onethird;
    rat_vector m_oneseventh;
};
#pragma GCC diagnostic pop

#endif /* UNLOGICALTESTCASE_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
