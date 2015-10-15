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

#define INFINT_STATICS_DECLARED 1

#include "exprtest_infint.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( ExprTestInfInt );

using namespace Commons::Math;

ExprTestInfInt::ExprTestInfInt() : CppUnit::TestFixture() {}

void ExprTestInfInt::setUp() {}

void ExprTestInfInt::tearDown() {}

void ExprTestInfInt::testExpression() {

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &a
#else
    const RationalExpressionTraits<infint_rational>::expr_type &a
#endif
    ( mk_rat_lit ( infint_rational ( 1, 8 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &b
#else
    const RationalExpressionTraits<infint_rational>::expr_type &b
#endif
    ( mk_rat_lit ( infint_rational ( 2, 73 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &c
#else
    const RationalExpressionTraits<infint_rational>::expr_type &c
#endif
    ( mk_rat_lit ( infint_rational ( 8, 17 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &d
#else
    const RationalExpressionTraits<infint_rational>::expr_type &d
#endif
    ( mk_rat_lit ( infint_rational ( 876, 2127 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &e
#else
    const RationalExpressionTraits<infint_rational>::expr_type &e
#endif
    ( mk_rat_lit ( infint_rational ( 670059l, 1656224l ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &f
#else
    const RationalExpressionTraits<infint_rational>::expr_type &f
#endif
    ( mk_rat_lit ( infint_rational ( -3, 2, 3 ) ) );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &g
#else
    const RationalExpressionTraits<infint_rational>::expr_type &g
#endif
    ( mk_rat_lit ( infint_rational ( 50.0 ) ) );

    const infint_rational &r_mod ( eval_rat_expr ( a % b ) );

    CPPUNIT_ASSERT_EQUAL ( 9l, r_mod.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 584l, r_mod.denominator().toLong() );

    const infint_rational &r_cpx ( eval_rat_expr ( ( ( ( ( a * b ) / -c ) % d ) - ( +e ) ) + f ) );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-7" ), r_cpx.numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 3l , r_cpx.denominator().toLong() );

    const infint_rational &r_cpx2 ( eval_rat_expr ( inv ( a * b / -c % d - e + f ) ) );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-3" ), r_cpx2.numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 7l, r_cpx2.denominator().toLong() );

    const infint_rational &r_cpx3 ( eval_rat_expr ( ( a * b ) / -c ) );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-17" ), r_cpx3.numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 2336l, r_cpx3.denominator().toLong() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &x ( mk_rat_proto_var ( g ) );
#else
    const RationalExpressionTraits<infint_rational>::variable_type &x ( mk_rat_proto_var ( g ) );
#endif

    const infint_rational &r0 ( eval_rat_expr ( x + g + x, infint_rational ( 2.0 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 54l, r0.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 1l, r0.denominator().toLong() );

    const infint_rational &r1 ( eval_rat_expr ( mk_rat_lit ( infint_rational ( 0, 1, 2 ) ) +
                                mk_rat_lit ( infint_rational ( 2, 3 ) ) +
                                mk_rat_lit ( infint_rational ( 3, 4 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 23l, r1.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r1.denominator().toLong() );

    const infint_rational &r2 ( eval_rat_expr ( mk_rat_lit ( r1 ) - infint_rational ( 23, 12 ) ) );

    CPPUNIT_ASSERT_EQUAL ( 0l, r2.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 1l, r2.denominator().toLong() );

    const infint_rational &r3 ( eval_rat_expr ( infint_rational ( 23, 12 ) -
                                mk_rat_lit ( infint_rational ( 22, 12 ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r3.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r3.denominator().toLong() );

    const infint_rational &r4 ( eval_rat_expr ( mk_rat_lit ( infint_rational ( 23, 12 ) ) +
                                ( -mk_rat_lit ( infint_rational ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r4.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r4.denominator().toLong() );

    const infint_rational &r5 ( eval_rat_expr ( infint_rational ( 23, 12 ) +
                                ( -mk_rat_lit ( infint_rational ( 22, 12 ) ) ) ) );

    CPPUNIT_ASSERT_EQUAL ( 1l, r5.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 12l, r5.denominator().toLong() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &expr ( a + b );
#endif

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const infint_rational &r6 ( eval_rat_expr ( expr + 1ul ) );
#else
    const infint_rational &r6 ( eval_rat_expr ( ( a + b ) + 1ul ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 673l , r6.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 584l , r6.denominator().toLong() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const infint_rational &r7 ( eval_rat_expr ( expr - 0.5 ) );
#else
    const infint_rational &r7 ( eval_rat_expr ( ( a + b ) - 0.5 ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-203" ), r7.numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 584l, r7.denominator().toLong() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const infint_rational &r8 ( eval_rat_expr ( expr * 0.5f ) );
#else
    const infint_rational &r8 ( eval_rat_expr ( ( a + b ) * 0.5f ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 89l, r8.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 1168l, r8.denominator().toLong() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const infint_rational &r9 ( eval_rat_expr ( expr / 8l ) );
#else
    const infint_rational &r9 ( eval_rat_expr ( ( a + b ) / 8l ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( 89l, r9.numerator().toLong() );
    CPPUNIT_ASSERT_EQUAL ( 4672l, r9.denominator().toLong() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const infint_rational &r10 ( eval_rat_expr ( expr % -0.3 ) );
#else
    const infint_rational &r10 ( eval_rat_expr ( ( a + b ) % -0.3 ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-431" ), r10.numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 2920l, r10.denominator().toLong() );

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const infint_rational &r11 ( eval_rat_expr ( inv ( expr % -0.3 ) ) );
#else
    const infint_rational &r11 ( eval_rat_expr ( inv ( ( a + b ) % -0.3 ) ) );
#endif

    CPPUNIT_ASSERT_EQUAL ( std::string ( "-2920" ), r11.numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( 431l, r11.denominator().toLong() );

}

void ExprTestInfInt::testIntegrate() {

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    const auto &x
#else
    const RationalExpressionTraits<infint_rational>::variable_type &x
#endif
    ( mk_rat_proto_var ( infint_rational() ) );

    const infint_rational &r ( integrate ( x / ( x + 1 ), 1, 5, 10 ) );

    CPPUNIT_ASSERT_EQUAL ( std::string ( "422563503196" ), r.numerator().toString() );
    CPPUNIT_ASSERT_EQUAL ( std::string ( "145568097675" ), r.denominator().toString() );

}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
