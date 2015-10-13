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

#ifndef CLNTESTCASE_H
#define CLNTESTCASE_H

#include <cppunit/extensions/HelperMacros.h>

#include "cln_rational.h"

CPPUNIT_NS_BEGIN

template<> struct assertion_traits<cln::cl_I> {

    static bool equal ( const cln::cl_I& x, const cln::cl_I& y ) {
        return x == y;
    }

    static std::string toString ( const cln::cl_I& x ) {
        std::ostringstream os;
        os << x;
        return os.str();
    }
};

CPPUNIT_NS_END

#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
class CLNTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( CLNTest );
    CPPUNIT_TEST ( testConstruct );
    CPPUNIT_TEST ( testConstructFromDouble );
    CPPUNIT_TEST ( testConstructFrom_cl_F_class );
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
    CPPUNIT_TEST_SUITE_END();

public:
    typedef cln::cl_I rational_type;

    CLNTest();

    void setUp();
    void tearDown();

    void testConstruct();
    void testConstructFromDouble();
    void testConstructFrom_cl_F_class();
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

private:
    typedef Commons::Math::cln_rational unchecked_sqrt;
    typedef std::vector<unchecked_sqrt> rat_vector_sqrt;
    typedef std::vector<Commons::Math::cln_rational> rat_vector;

    unchecked_sqrt m_sqrt2;
    rat_vector_sqrt m_twosqrt;
    rat_vector m_onethird;
    rat_vector m_oneseventh;
};
#pragma GCC diagnostic pop

#endif /* CLNTESTCASE_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
