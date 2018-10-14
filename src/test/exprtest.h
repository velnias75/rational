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

#ifndef EXPRTESTCASE_H
#define EXPRTESTCASE_H

#if defined(HAVE_CONFIG_H) || defined(IN_IDE_PARSER)
#include "config.h"
#endif

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
#define RATIONAL_OVERRIDE override
#define RATIONAL_FINAL final
#else
#define RATIONAL_OVERRIDE
#define RATIONAL_FINAL
#endif

#include <cppunit/extensions/HelperMacros.h>

#include "expr_rational.h"

#ifdef HAVE_GMPXX_H
#include "gmp_rational.h"
#endif

#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
class ExprTest RATIONAL_FINAL : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( ExprTest );
    CPPUNIT_TEST ( testExpression );
    CPPUNIT_TEST ( testExpression_gmp );
    CPPUNIT_TEST ( testIntegrate );
    CPPUNIT_TEST_SUITE_END();

public:
    typedef Commons::Math::Rational<long, Commons::Math::GCD_euclid_fast,
            Commons::Math::NO_OPERATOR_CHECK> long_rational;

    ExprTest();

    void setUp() RATIONAL_OVERRIDE;
    void tearDown() RATIONAL_OVERRIDE;

    void testExpression();
    void testExpression_gmp();
    void testIntegrate();

private:
#ifdef HAVE_GMPXX_H
    template<class ExprT> inline static
    Commons::Math::gmp_rational integrate ( const ExprT &e, const Commons::Math::gmp_rational &from,
                                            const Commons::Math::gmp_rational &to,
                                            std::size_t n ) {
        Commons::Math::gmp_rational sum;
        const static Commons::Math::gmp_rational two ( 2, 1 );
        const Commons::Math::gmp_rational &step ( ( to - from ) / n );

        for ( Commons::Math::gmp_rational i ( from + ( step / two ) ); i < to; i += step ) {
            sum += Commons::Math::eval_rat_expr ( e, i );
        }

        return step * sum;
    }
#endif

};
#pragma GCC diagnostic pop

#endif /* EXPRTESTCASE_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
