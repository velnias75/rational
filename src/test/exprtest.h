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

#ifndef EXPRTESTCASE_H
#define EXPRTESTCASE_H

#if defined(HAVE_CONFIG_H) || defined(IN_IDE_PARSER)
#include "config.h"
#endif

#include <cppunit/extensions/HelperMacros.h>

#include "expr_rational.h"

#ifdef HAVE_GMPXX_H
#include "gmp_rational.h"
#endif

#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
class ExprTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( ExprTest );
    CPPUNIT_TEST ( testExpression );
    CPPUNIT_TEST ( testExpression_gmp );
    CPPUNIT_TEST_SUITE_END();

public:
    ExprTest();

    void setUp();
    void tearDown();

    void testExpression();
    void testExpression_gmp();
};
#pragma GCC diagnostic pop

#endif /* EXPRTESTCASE_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
