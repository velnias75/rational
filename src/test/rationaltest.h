/*
 * Copyright 2015-2018 by Heiko Schäfer <heiko@rangun.de>
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

#ifndef RATIONALTESTCASE_H
#define RATIONALTESTCASE_H

#include <vector>

#include <stdint.h>

#include <cppunit/extensions/HelperMacros.h>

#include "rational.h"

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
class RationalTest RATIONAL_FINAL : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( RationalTest );
    CPPUNIT_TEST ( testNullRational );
    CPPUNIT_TEST ( testConstruct );
    CPPUNIT_TEST ( testConstructFromDouble );
    CPPUNIT_TEST ( testConstructFromExpression );
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
    CPPUNIT_TEST ( testString );
    CPPUNIT_TEST ( testIOStreamOps );
    CPPUNIT_TEST ( testPrecision );
    CPPUNIT_TEST ( testAlgorithm );
    CPPUNIT_TEST ( testStdMath );
    CPPUNIT_TEST ( testRatRat );
    CPPUNIT_TEST ( testGoldenRatio );
    CPPUNIT_TEST_SUITE_END();

public:
    typedef int32_t rational_type;

    RationalTest();

    void setUp() RATIONAL_OVERRIDE;
    void tearDown() RATIONAL_OVERRIDE;

    void testNullRational();
    void testConstruct();
    void testConstructFromDouble();
    void testConstructFromExpression();
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
    void testString();
    void testIOStreamOps();
    void testPrecision();
    void testAlgorithm();
    void testStdMath();
    void testRatRat();
    void testGoldenRatio();

private:
    typedef Commons::Math::Rational<uint64_t, Commons::Math::GCD_euclid,
            Commons::Math::ENABLE_OVERFLOW_CHECK> checked_sqrt;

    Commons::Math::Rational<rational_type> m_nullRational;
    checked_sqrt m_sqrt2;

    typedef std::vector<Commons::Math::Rational<rational_type> > rat_vector;
    typedef std::vector<Commons::Math::Rational<rational_type,
            Commons::Math::GCD_stein> > rat_vector_stein;

    rat_vector m_accu;
    rat_vector_stein m_accu_stein;
    rat_vector m_onethird;
    rat_vector m_oneseventh;

    typedef std::vector<Commons::Math::Rational<uint64_t> > rat_vector_ul;
    rat_vector_ul m_accu_ul;

    typedef std::vector<checked_sqrt> rat_vector_sqrt;
    rat_vector_sqrt m_twosqrt;
};
#pragma GCC diagnostic pop

#endif /* RATIONALTESTCASE_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
