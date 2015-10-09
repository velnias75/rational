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

    const RationalExpressionTraits<Rational<long> >::expr_type &a (
        mk_rat_lit ( Rational<long> ( 1, 8 ) ) );

    const RationalExpressionTraits<Rational<long> >::expr_type &b (
        mk_rat_lit ( Rational<long> ( 2, 73 ) ) );

    const RationalExpressionTraits<Rational<long> >::expr_type &c (
        mk_rat_lit ( Rational<long> ( 8, 17 ) ) );

    const RationalExpressionTraits<Rational<long> >::expr_type &d (
        mk_rat_lit ( Rational<long> ( 876, 2127 ) ) );

    const RationalExpressionTraits<Rational<long> >::expr_type &e (
        mk_rat_lit ( Rational<long> ( 670059l, 1656224l ) ) );

    const RationalExpressionTraits<Rational<long> >::expr_type &f (
        mk_rat_lit ( Rational<long> ( -3, 2, 3 ) ) );

    const Rational<long> r_mod ( eval_rat_expr ( a % b ) );

    CPPUNIT_ASSERT_EQUAL ( 9l, r_mod.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 584l, r_mod.denominator() );

    const Rational<long> r_cpx ( eval_rat_expr ( ( ( ( ( a * b ) / -c ) % d ) - ( +e ) ) + f ) );

    CPPUNIT_ASSERT_EQUAL ( -7l, r_cpx.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 3l, r_cpx.denominator() );

    const Rational<long> r_cpx2 ( eval_rat_expr ( inv ( a * b / -c % d - e + f ) ) );

    CPPUNIT_ASSERT_EQUAL ( -3l, r_cpx2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 7l, r_cpx2.denominator() );

    const Rational<long> r_cpx3 ( eval_rat_expr ( ( ( a * b ) / -c ) ) );

    CPPUNIT_ASSERT_EQUAL ( -17l, r_cpx3.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 2336l, r_cpx3.denominator() );

    const RationalExpressionTraits<Rational<long> >::expr_type &l (
        mk_rat_lit ( Rational<long> ( 50.0 ) ) );

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

    const Rational<long> r4 ( eval_rat_expr ( mk_rat_lit ( Rational<long> ( 23, 12 ) ) +
                              ( -mk_rat_lit ( Rational<long> ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r4.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r4.denominator() );

    const Rational<long> r5 ( eval_rat_expr ( Rational<long> ( 23, 12 ) +
                              ( -mk_rat_lit ( Rational<long> ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r5.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r5.denominator() );
}

void ExprTest::testExpression_gmp() {
#ifdef HAVE_GMPXX_H

    const RationalExpressionTraits<gmp_rational>::expr_type &a (
        mk_rat_lit ( gmp_rational ( 1, 8 ) ) );

    const RationalExpressionTraits<gmp_rational>::expr_type &b (
        mk_rat_lit ( gmp_rational ( 2, 73 ) ) );

    const RationalExpressionTraits<gmp_rational>::expr_type &c (
        mk_rat_lit ( gmp_rational ( 8, 17 ) ) );

    const RationalExpressionTraits<gmp_rational>::expr_type &d (
        mk_rat_lit ( gmp_rational ( 876, 2127 ) ) );

    const RationalExpressionTraits<gmp_rational>::expr_type &e (
        mk_rat_lit ( gmp_rational ( 670059l, 1656224l ) ) );

    const RationalExpressionTraits<gmp_rational>::expr_type &f (
        mk_rat_lit ( gmp_rational ( -3, 2, 3 ) ) );

    const gmp_rational r_mod ( eval_rat_expr ( a % b ) );

    CPPUNIT_ASSERT_EQUAL ( 9l, r_mod.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 584l, r_mod.denominator().get_si() );

    const gmp_rational r_cpx ( eval_rat_expr ( ( ( ( ( a * b ) / -c ) % d ) - ( +e ) ) + f ) );

    CPPUNIT_ASSERT_EQUAL ( -7l, r_cpx.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 3l, r_cpx.denominator().get_si() );

    const gmp_rational r_cpx2 ( eval_rat_expr ( a * b / -c % d - e + f ) );

    CPPUNIT_ASSERT_EQUAL ( -7l, r_cpx2.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 3l, r_cpx2.denominator().get_si() );

    const gmp_rational r_cpx3 ( eval_rat_expr ( ( a * b ) / -c ) );

    CPPUNIT_ASSERT_EQUAL ( -17l, r_cpx3.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 2336l, r_cpx3.denominator().get_si() );

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

    const gmp_rational r4 ( eval_rat_expr ( mk_rat_lit ( gmp_rational ( 23, 12 ) ) +
                                            ( -mk_rat_lit ( gmp_rational ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r4.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r4.denominator().get_si() );

    const gmp_rational r5 ( eval_rat_expr ( gmp_rational ( 23, 12 ) +
                                            ( -mk_rat_lit ( gmp_rational ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r5.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r5.denominator().get_si() );

#endif
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
