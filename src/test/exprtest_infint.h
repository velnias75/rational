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

#ifndef EXPRTESTCASE_INFINT_H
#define EXPRTESTCASE_INFINT_H

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

#include "infint_rational.h"

#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
class ExprTestInfInt RATIONAL_FINAL : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( ExprTestInfInt );
    CPPUNIT_TEST ( testExpression );
    CPPUNIT_TEST ( testIntegrate );
    CPPUNIT_TEST_SUITE_END();

public:
    ExprTestInfInt();

    void setUp() RATIONAL_OVERRIDE;
    void tearDown() RATIONAL_OVERRIDE;

    void testExpression();
    void testIntegrate();

private:
    template<class ExprT> inline static
    Commons::Math::infint_rational integrate ( const ExprT &e,
            const Commons::Math::infint_rational &from,
            const Commons::Math::infint_rational &to, std::size_t n ) {

        Commons::Math::infint_rational sum;
        const static Commons::Math::infint_rational two ( 2, 1 );
        const Commons::Math::infint_rational &step ( ( to - from ) / n );

        for ( Commons::Math::infint_rational i ( from + ( step / two ) ); i < to; i += step ) {
            sum += Commons::Math::eval_rat_expr ( e, i );
        }

        return step * sum;
    }

};
#pragma GCC diagnostic pop

#endif /* EXPRTESTCASE_INFINT_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
