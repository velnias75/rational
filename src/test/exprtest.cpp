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

    declare_rat_var ( long, x );

    RationalExpressionTraits<Rational<long> >::expr_type l ( mk_rat_lit
            ( Rational<long> ( 50.0 ) ) );

    const Rational<long> r1 ( eval_rat_exp ( x + l + x, Rational<long> ( 2.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 54l, r1.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1l, r1.denominator() );

    const Rational<long> r2 ( eval_rat_exp ( mk_rat_lit ( Rational<long> ( 1.0/2.0 ) ) +
                              mk_rat_lit ( Rational<long> ( 2.0/3.0 ) ) +
                              mk_rat_lit ( Rational<long> ( 3.0/4.0 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 23l, r2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r2.denominator() );

    const Rational<long> r3 ( eval_rat_exp ( Rational<long> ( 1.0/2.0 ) +
                              Rational<long> ( 2.0/3.0 ) + Rational<long> ( 3.0/4.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 23l, r3.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r3.denominator() );

    const gmp_rational r4 ( eval_rat_exp ( mk_rat_lit ( gmp_rational ( 0, 1, 2 ) ) +
                                           mk_rat_lit ( gmp_rational ( 2, 3 ) ) +
                                           mk_rat_lit ( gmp_rational ( 3, 4 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 23l, r4.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r4.denominator().get_si() );

    const gmp_rational r5 ( eval_rat_exp ( mk_rat_lit ( r4 ) - gmp_rational ( 23, 12 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 0l, r5.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, r5.denominator().get_si() );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
