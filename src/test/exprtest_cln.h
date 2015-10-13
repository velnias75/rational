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

#ifndef EXPRTESTCASE_CLN_H
#define EXPRTESTCASE_CLN_H

#if defined(HAVE_CONFIG_H) || defined(IN_IDE_PARSER)
#include "config.h"
#endif

#include <cppunit/extensions/HelperMacros.h>

#include "expr_rational.h"

#include "cln_rational.h"

#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
class ExprTestCLN : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( ExprTestCLN );
    CPPUNIT_TEST ( testExpression );
    CPPUNIT_TEST ( testIntegrate );
    CPPUNIT_TEST_SUITE_END();

public:
    ExprTestCLN();

    void setUp();
    void tearDown();

    void testExpression();
    void testIntegrate();

private:
    template<class ExprT> inline static
    Commons::Math::cln_rational integrate ( const ExprT &e, const Commons::Math::cln_rational &from,
                                            const Commons::Math::cln_rational &to,
                                            std::size_t n ) {
        Commons::Math::cln_rational sum;
        const static Commons::Math::cln_rational two ( 2, 1 );
        const Commons::Math::cln_rational &step ( ( to - from ) / n );

        for ( Commons::Math::cln_rational i ( from + ( step / two ) ); i < to; i += step ) {
            sum += Commons::Math::eval_rat_expr ( e, i );
        }

        return step * sum;
    }

};
#pragma GCC diagnostic pop

#endif /* EXPRTESTCASE_CLN_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 