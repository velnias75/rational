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

#include "exprtest.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( ExprTest );

using namespace Commons::Math;

ExprTest::ExprTest() : CppUnit::TestFixture() {}

void ExprTest::setUp() {}

void ExprTest::tearDown() {}

void ExprTest::testExpression() {

    CPPUNIT_ASSERT_EQUAL ( 23l, ( mk_rat_lit<gmp_rational::integer_type> ( 0, 1, 2 ) +
                                  mk_rat_lit<gmp_rational::integer_type> ( 2, 3 ) +
                                  mk_rat_lit<gmp_rational::integer_type> ( 3, 4 ) ) ()
                           .numerator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 12l, ( mk_rat_lit ( gmp_rational ( 1, 2 ) ) +
                                  mk_rat_lit ( gmp_rational ( 2, 3 ) ) +
                                  mk_rat_lit ( gmp_rational ( 3, 4 ) ) ) ()
                           .denominator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 0l, ( ( mk_rat_lit<gmp_rational::integer_type> ( 0, 1, 2 ) +
                                   mk_rat_lit<gmp_rational::integer_type> ( 2, 3 ) +
                                   mk_rat_lit<gmp_rational::integer_type> ( 3, 4 ) ) -
                                 mk_rat_lit<gmp_rational::integer_type> ( 23, 12 ) ) ()
                           .numerator().get_si() );

    CPPUNIT_ASSERT_EQUAL ( 1l, ( ( mk_rat_lit ( gmp_rational ( 1, 2 ) ) +
                                   mk_rat_lit ( gmp_rational ( 2, 3 ) ) +
                                   mk_rat_lit ( gmp_rational ( 3, 4 ) ) ) -
                                 mk_rat_lit ( gmp_rational ( 23, 12 ) ) ) ()
                           .denominator().get_si() );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
