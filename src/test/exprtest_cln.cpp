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

#include "exprtest_cln.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( ExprTestCLN );

using namespace Commons::Math;

ExprTestCLN::ExprTestCLN() : CppUnit::TestFixture() {}

void ExprTestCLN::setUp() {}

void ExprTestCLN::tearDown() {}

void ExprTestCLN::testExpression() {

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &a
#else
    const RationalExpressionTraits<cln_rational>::expr_type &a
#endif
    ( mk_rat_lit ( cln_rational ( 1, 8 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &b
#else
    const RationalExpressionTraits<cln_rational>::expr_type &b
#endif
    ( mk_rat_lit ( cln_rational ( 2, 73 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &c
#else
    const RationalExpressionTraits<cln_rational>::expr_type &c
#endif
    ( mk_rat_lit ( cln_rational ( 8, 17 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &d
#else
    const RationalExpressionTraits<cln_rational>::expr_type &d
#endif
    ( mk_rat_lit ( cln_rational ( 876, 2127 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &e
#else
    const RationalExpressionTraits<cln_rational>::expr_type &e
#endif
    ( mk_rat_lit ( cln_rational ( 670059l, 1656224l ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &f
#else
    const RationalExpressionTraits<cln_rational>::expr_type &f
#endif
    ( mk_rat_lit ( cln_rational ( -3, 2, 3 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &g
#else
    const RationalExpressionTraits<cln_rational>::expr_type &g
#endif
    ( mk_rat_lit ( cln_rational ( 50.0 ) ) );

    const cln_rational &r_mod ( eval_rat_expr ( a % b ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 9l ), r_mod.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 584l ), r_mod.denominator() );

    const cln_rational &r_cpx ( eval_rat_expr ( ( ( ( ( a * b ) / -c ) % d ) - ( +e ) ) + f ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -7l ), r_cpx.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 3l ), r_cpx.denominator() );

    const cln_rational &r_cpx2 ( eval_rat_expr ( inv ( a * b / -c % d - e + f ) ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -3l ), r_cpx2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 7l ), r_cpx2.denominator() );

    const cln_rational &r_cpx3 ( eval_rat_expr ( ( a * b ) / -c ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -17l ), r_cpx3.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2336l ), r_cpx3.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &x ( mk_rat_proto_var ( g ) );
#else
    const RationalExpressionTraits<cln_rational>::variable_type &x ( mk_rat_proto_var ( g ) );
#endif

    const cln_rational &r0 ( eval_rat_expr ( x + g + x, cln_rational ( 2.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 54l ), r0.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), r0.denominator() );

    const cln_rational &r1 ( eval_rat_expr ( mk_rat_lit ( cln_rational ( 0, 1, 2 ) ) +
                             mk_rat_lit ( cln_rational ( 2, 3 ) ) +
                             mk_rat_lit ( cln_rational ( 3, 4 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 23l ), r1.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 12l ), r1.denominator() );

    const cln_rational &r2 ( eval_rat_expr ( mk_rat_lit ( r1 ) - cln_rational ( 23, 12 ) ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 0l ), r2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), r2.denominator() );

    const cln_rational &r3 ( eval_rat_expr ( cln_rational ( 23, 12 ) -
                             mk_rat_lit ( cln_rational ( 22, 12 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), r3.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 12l ), r3.denominator() );

    const cln_rational &r4 ( eval_rat_expr ( mk_rat_lit ( cln_rational ( 23, 12 ) ) +
                             ( -mk_rat_lit ( cln_rational ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), r4.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 12l ), r4.denominator() );

    const cln_rational &r5 ( eval_rat_expr ( cln_rational ( 23, 12 ) +
                             ( -mk_rat_lit ( cln_rational ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1l ), r5.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 12l ), r5.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &expr ( a + b );
#endif

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const cln_rational &r6 ( eval_rat_expr ( expr + 1ul ) );
#else
    const cln_rational &r6 ( eval_rat_expr ( ( a + b ) + 1ul ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 673l ), r6.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 584l ), r6.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const cln_rational &r7 ( eval_rat_expr ( expr - 0.5 ) );
#else
    const cln_rational &r7 ( eval_rat_expr ( ( a + b ) - 0.5 ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -203l ), r7.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 584l ), r7.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const cln_rational &r8 ( eval_rat_expr ( expr * 0.5f ) );
#else
    const cln_rational &r8 ( eval_rat_expr ( ( a + b ) * 0.5f ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 89l ), r8.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 1168l ), r8.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const cln_rational &r9 ( eval_rat_expr ( expr / 8l ) );
#else
    const cln_rational &r9 ( eval_rat_expr ( ( a + b ) / 8l ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 89l ), r9.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 4672l ), r9.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const cln_rational &r10 ( eval_rat_expr ( expr % -0.3 ) );
#else
    const cln_rational &r10 ( eval_rat_expr ( ( a + b ) % -0.3 ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -431l ), r10.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 2920l ), r10.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const cln_rational &r11 ( eval_rat_expr ( inv ( expr % -0.3 ) ) );
#else
    const cln_rational &r11 ( eval_rat_expr ( inv ( ( a + b ) % -0.3 ) ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( -2920l ), r11.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( 431l ), r11.denominator() );

}

void ExprTestCLN::testIntegrate() {

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &x
#else
    const RationalExpressionTraits<cln_rational>::variable_type &x
#endif
    ( mk_rat_proto_var ( cln_rational() ) );

    const cln_rational &r ( integrate ( x / ( 1 + x ), 1, 5, 10 ) );

    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "422563503196" ), r.numerator() );
    CPPUNIT_ASSERT_EQUAL ( cln::cl_I ( "145568097675" ), r.denominator() );

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
