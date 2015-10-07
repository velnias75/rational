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

    const RationalExpressionTraits<Rational<long> >::expr_type &l ( mk_rat_lit
            ( Rational<long> ( 50.0 ) ) );

    const RationalExpressionTraits<Rational<long> >::variable_type &x ( mk_rat_proto_var ( l ) );

    const Rational<long> r1 ( eval_rat_expr ( x + l + x, Rational<long> ( 2.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 54l, r1.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1l, r1.denominator() );

    const Rational<long> r2 ( eval_rat_expr ( mk_rat_lit ( Rational<long> ( 1.0/2.0 ) ) +
                              mk_rat_lit ( Rational<long> ( 2.0/3.0 ) ) +
                              mk_rat_lit ( Rational<long> ( 3.0/4.0 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 23l, r2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r2.denominator() );

    const Rational<long> r3 ( eval_rat_expr ( Rational<long> ( 1.0/2.0 ) +
                              Rational<long> ( 2.0/3.0 ) + Rational<long> ( 3.0/4.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 23l, r3.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r3.denominator() );
}

void ExprTest::testExpression_gmp() {
#ifdef HAVE_GMPXX_H

    const gmp_rational r1 ( eval_rat_expr ( mk_rat_lit ( gmp_rational ( 0, 1, 2 ) ) +
                                            mk_rat_lit ( gmp_rational ( 2, 3 ) ) +
                                            mk_rat_lit ( gmp_rational ( 3, 4 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 23l, r1.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r1.denominator().get_si() );

    const gmp_rational r2 ( eval_rat_expr ( mk_rat_lit ( r1 ) - gmp_rational ( 23, 12 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 0l, r2.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, r2.denominator().get_si() );

    const gmp_rational r3 ( eval_rat_expr ( gmp_rational ( 23, 12 ) -
                                            mk_rat_lit ( gmp_rational ( 22, 12 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r3.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r3.denominator().get_si() );

#endif
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
