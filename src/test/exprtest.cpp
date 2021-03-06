/*
 * Copyright 2015-2016 by Heiko Schäfer <heiko@rangun.de>
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

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &a
#else
    const RationalExpressionTraits<long_rational>::expr_type &a
#endif
    ( mk_rat_lit ( long_rational ( 1, 8 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &b
#else
    const RationalExpressionTraits<long_rational>::expr_type &b
#endif
    ( mk_rat_lit ( long_rational ( 2, 73 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &c
#else
    const RationalExpressionTraits<long_rational>::expr_type &c
#endif
    ( mk_rat_lit ( long_rational ( 8, 17 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &d
#else
    const RationalExpressionTraits<long_rational>::expr_type &d
#endif
    ( mk_rat_lit ( long_rational ( 876, 2127 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &e
#else
    const RationalExpressionTraits<long_rational>::expr_type &e
#endif
    ( mk_rat_lit ( long_rational ( 670059l, 1656224l ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &f
#else
    const RationalExpressionTraits<long_rational>::expr_type &f
#endif
    ( mk_rat_lit ( long_rational ( -3, 2, 3 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &g
#else
    const RationalExpressionTraits<long_rational>::expr_type &g
#endif
    ( mk_rat_lit ( long_rational ( 50.0 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &h
#else
    const RationalExpressionTraits<long_rational>::expr_type &h
#endif
    ( mk_rat_lit ( long_rational ( 85, -48 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &i
#else
    const RationalExpressionTraits<long_rational>::expr_type &i
#endif
    ( mk_rat_lit ( long_rational ( 7, 3 ) ) );

    const long_rational &r_mod ( eval_rat_expr ( a % b ) );

    CPPUNIT_ASSERT_EQUAL ( 9l, r_mod.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 584l, r_mod.denominator() );

    const long_rational &r_cpx ( eval_rat_expr ( ( ( ( ( a * b ) / -c ) % d ) - ( +e ) ) + f ) );

    CPPUNIT_ASSERT_EQUAL ( -7l, r_cpx.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 3l, r_cpx.denominator() );

    const long_rational &r_cpx2 ( eval_rat_expr ( inv ( a * b / -c % d - e + f ) ) );

    CPPUNIT_ASSERT_EQUAL ( -3l, r_cpx2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 7l, r_cpx2.denominator() );

    const long_rational &r_cpx3 ( eval_rat_expr ( ( ( a * b ) / -c ) ) );

    CPPUNIT_ASSERT_EQUAL ( -17l, r_cpx3.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 2336l, r_cpx3.denominator() );

    const RationalExpressionTraits<long_rational>::variable_type &x ( mk_rat_proto_var ( g ) );

    const long_rational &r1 ( eval_rat_expr ( x + g + x, long_rational ( 2.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 54l, r1.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1l, r1.denominator() );

    const long_rational &r2 ( eval_rat_expr ( mk_rat_lit ( long_rational ( 1.0/2.0 ) ) +
                              mk_rat_lit ( long_rational ( 2.0/3.0 ) ) +
                              mk_rat_lit ( long_rational ( 3.0/4.0 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 23l, r2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r2.denominator() );

    const long_rational &r3 ( eval_rat_expr ( long_rational ( 1.0/2.0 ) +
                              long_rational ( 2.0/3.0 ) + long_rational ( 3.0/4.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 23l, r3.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r3.denominator() );

    const long_rational &r4 ( eval_rat_expr ( mk_rat_lit ( long_rational ( 23, 12 ) ) +
                              ( -mk_rat_lit ( long_rational ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r4.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r4.denominator() );

    const long_rational &r5 ( eval_rat_expr ( long_rational ( 23, 12 ) +
                              ( -mk_rat_lit ( long_rational ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r5.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r5.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &expr ( a + b );
#endif

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const long_rational &r6 ( eval_rat_expr ( 1.0f + expr + 1ul ) );
#else
    const long_rational &r6 ( eval_rat_expr ( 1.0f + ( a + b ) + 1ul ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 1257l, r6.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 584l, r6.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const long_rational &r61 ( eval_rat_expr ( 1.0f - expr + 1ul ) );
#else
    const long_rational &r61 ( eval_rat_expr ( 1.0f - ( a + b ) + 1ul ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 1079l, r61.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 584l, r61.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const long_rational &r62 ( eval_rat_expr ( 1.0 / expr + 1ul ) );
#else
    const long_rational &r62 ( eval_rat_expr ( 1.0 / ( a + b ) + 1ul ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 673l, r62.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 89l, r62.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const long_rational &r63 ( eval_rat_expr ( 8 * expr + 1ul ) );
#else
    const long_rational &r63 ( eval_rat_expr ( 8 * ( a + b ) + 1ul ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 162l, r63.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 73l, r63.denominator() );

    CPPUNIT_ASSERT_EQUAL ( 673l, r62.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 89l, r62.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const long_rational &r64 ( eval_rat_expr ( 8u % ( expr + 1ul ) ) );
#else
    const long_rational &r64 ( eval_rat_expr ( 8u % ( ( a + b ) + 1ul ) ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 317l, r64.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 292l, r64.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const long_rational &r7 ( eval_rat_expr ( expr - 0.5 ) );
#else
    const long_rational &r7 ( eval_rat_expr ( ( a + b ) - 0.5 ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( -203l, r7.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 584l, r7.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const long_rational &r8 ( eval_rat_expr ( expr * 0.5f ) );
#else
    const long_rational &r8 ( eval_rat_expr ( ( a + b ) * 0.5f ) );
#endif
    CPPUNIT_ASSERT_EQUAL ( 89l, r8.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 1168l, r8.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const long_rational &r9 ( eval_rat_expr ( expr / 8l ) );
#else
    const long_rational &r9 ( eval_rat_expr ( ( a + b ) / 8l ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 89l, r9.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4672l, r9.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const long_rational &r10 ( eval_rat_expr ( expr % -0.3l ) );
#else
    const long_rational &r10 ( eval_rat_expr ( ( a + b ) % -0.3l ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( -431l, r10.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 2920l, r10.denominator() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const long_rational &r11 ( eval_rat_expr ( inv ( expr % -0.3l ) ) );
#else
    const long_rational &r11 ( eval_rat_expr ( inv ( ( a + b ) % -0.3l ) ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( -2920l, r11.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 431l, r11.denominator() );

    const long_rational &sqr1 ( eval_rat_expr ( sqrt ( h + i ) ) );

    CPPUNIT_ASSERT_EQUAL ( 3l, sqr1.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4l, sqr1.denominator() );

    const long_rational &sqr2 ( eval_rat_expr ( sqrt ( i + h ) ) );

    CPPUNIT_ASSERT_EQUAL ( 3l, sqr2.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 4l, sqr2.denominator() );

    const long_rational &sqr3 ( eval_rat_expr ( sqrt ( h + i ) * sqrt ( i + h ) ) );

    CPPUNIT_ASSERT_EQUAL ( 9l, sqr3.numerator() );
    CPPUNIT_ASSERT_EQUAL ( 16l, sqr3.denominator() );
}

void ExprTest::testExpression_gmp() {
#ifdef HAVE_GMPXX_H

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &a
#else
    const RationalExpressionTraits<gmp_rational>::expr_type &a
#endif
    ( mk_rat_lit ( gmp_rational ( 1, 8 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &b
#else
    const RationalExpressionTraits<gmp_rational>::expr_type &b
#endif
    ( mk_rat_lit ( gmp_rational ( 2, 73 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &c
#else
    const RationalExpressionTraits<gmp_rational>::expr_type &c
#endif
    ( mk_rat_lit ( gmp_rational ( 8, 17 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &d
#else
    const RationalExpressionTraits<gmp_rational>::expr_type &d
#endif
    ( mk_rat_lit ( gmp_rational ( 876, 2127 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &e
#else
    const RationalExpressionTraits<gmp_rational>::expr_type &e
#endif
    ( mk_rat_lit ( gmp_rational ( 670059l, 1656224l ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &f
#else
    const RationalExpressionTraits<gmp_rational>::expr_type &f
#endif
    ( mk_rat_lit ( gmp_rational ( -3, 2, 3 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &g
#else
    const RationalExpressionTraits<gmp_rational>::expr_type &g
#endif
    ( mk_rat_lit ( gmp_rational ( 50.0 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &h
#else
    const RationalExpressionTraits<gmp_rational>::expr_type &h
#endif
    ( mk_rat_lit ( gmp_rational ( 85, -48 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &i
#else
    const RationalExpressionTraits<gmp_rational>::expr_type &i
#endif
    ( mk_rat_lit ( gmp_rational ( 7, 3 ) ) );

    const gmp_rational &r_mod ( eval_rat_expr ( a % b ) );

    CPPUNIT_ASSERT_EQUAL ( 9l, r_mod.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 584l, r_mod.denominator().get_si() );

    const gmp_rational &r_cpx ( eval_rat_expr ( ( ( ( ( a * b ) / -c ) % d ) - ( +e ) ) + f ) );

    CPPUNIT_ASSERT_EQUAL ( -7l, r_cpx.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 3l, r_cpx.denominator().get_si() );

    const gmp_rational &r_cpx2 ( eval_rat_expr ( inv ( a * b / -c % d - e + f ) ) );

    CPPUNIT_ASSERT_EQUAL ( -3l, r_cpx2.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 7l, r_cpx2.denominator().get_si() );

    const gmp_rational &r_cpx3 ( eval_rat_expr ( ( a * b ) / -c ) );

    CPPUNIT_ASSERT_EQUAL ( -17l, r_cpx3.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 2336l, r_cpx3.denominator().get_si() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &x ( mk_rat_proto_var ( g ) );
#else
    const RationalExpressionTraits<gmp_rational>::variable_type &x ( mk_rat_proto_var ( g ) );
#endif

    const gmp_rational &r0 ( eval_rat_expr ( x + g + x, gmp_rational ( 2.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 54l, r0.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, r0.denominator().get_si() );

    const gmp_rational &r1 ( eval_rat_expr ( mk_rat_lit ( gmp_rational ( 0, 1, 2 ) ) +
                             mk_rat_lit ( gmp_rational ( 2, 3 ) ) +
                             mk_rat_lit ( gmp_rational ( 3, 4 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 23l, r1.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r1.denominator().get_si() );

    const gmp_rational &r2 ( eval_rat_expr ( mk_rat_lit ( r1 ) - gmp_rational ( 23, 12 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 0l, r2.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1l, r2.denominator().get_si() );

    const gmp_rational &r3 ( eval_rat_expr ( gmp_rational ( 23, 12 ) -
                             mk_rat_lit ( gmp_rational ( 22, 12 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r3.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r3.denominator().get_si() );

    const gmp_rational &r4 ( eval_rat_expr ( mk_rat_lit ( gmp_rational ( 23, 12 ) ) +
                             ( -mk_rat_lit ( gmp_rational ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r4.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r4.denominator().get_si() );

    const gmp_rational &r5 ( eval_rat_expr ( gmp_rational ( 23, 12 ) +
                             ( -mk_rat_lit ( gmp_rational ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r5.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r5.denominator().get_si() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &expr ( a + b );
#endif

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const gmp_rational &r6 ( eval_rat_expr ( expr + 1ul ) );
#else
    const gmp_rational &r6 ( eval_rat_expr ( ( a + b ) + 1ul ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 673l, r6.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 584l, r6.denominator().get_si() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const gmp_rational &r7 ( eval_rat_expr ( expr - 0.5 ) );
#else
    const gmp_rational &r7 ( eval_rat_expr ( ( a + b ) - 0.5 ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( -203l, r7.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 584l, r7.denominator().get_si() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const gmp_rational &r8 ( eval_rat_expr ( expr * 0.5f ) );
#else
    const gmp_rational &r8 ( eval_rat_expr ( ( a + b ) * 0.5f ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 89l, r8.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 1168l, r8.denominator().get_si() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const gmp_rational &r9 ( eval_rat_expr ( expr / 8l ) );
#else
    const gmp_rational &r9 ( eval_rat_expr ( ( a + b ) / 8l ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 89l, r9.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4672l, r9.denominator().get_si() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const gmp_rational &r10 ( eval_rat_expr ( expr % -0.3 ) );
#else
    const gmp_rational &r10 ( eval_rat_expr ( ( a + b ) % -0.3 ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( -431l, r10.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 2920l, r10.denominator().get_si() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const gmp_rational &r11 ( eval_rat_expr ( inv ( expr % -0.3 ) ) );
#else
    const gmp_rational &r11 ( eval_rat_expr ( inv ( ( a + b ) % -0.3 ) ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( -2920l, r11.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 431l, r11.denominator().get_si() );

    const gmp_rational &sqr1 ( eval_rat_expr ( sqrt ( h + i ) ) );

    CPPUNIT_ASSERT_EQUAL ( 3l, sqr1.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4l, sqr1.denominator().get_si() );

    const gmp_rational &sqr2 ( eval_rat_expr ( sqrt ( i + h ) ) );

    CPPUNIT_ASSERT_EQUAL ( 3l, sqr2.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 4l, sqr2.denominator().get_si() );

    const gmp_rational &sqr3 ( eval_rat_expr ( sqrt ( h + i ) * sqrt ( i + h ) ) );

    CPPUNIT_ASSERT_EQUAL ( 9l, sqr3.numerator().get_si() );
    CPPUNIT_ASSERT_EQUAL ( 16l, sqr3.denominator().get_si() );
#endif
}

void ExprTest::testIntegrate() {
#ifdef HAVE_GMPXX_H

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &x
#else
    const RationalExpressionTraits<gmp_rational>::variable_type &x
#endif
    ( mk_rat_proto_var ( gmp_rational() ) );

    const gmp_rational &r ( integrate ( x / ( 1 + x ), 1, 5, 10 ) );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "422563503196" ), r.numerator().get_str() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "145568097675" ), r.denominator().get_str() );

#endif
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
