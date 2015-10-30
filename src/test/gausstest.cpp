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

#include "gausstest.h"

CPPUNIT_TEST_SUITE_REGISTRATION ( GaussTest );

using namespace Commons::Math;

GaussTest::GaussTest() : CppUnit::TestFixture() {}

void GaussTest::setUp() {}

void GaussTest::tearDown() {}

void GaussTest::testConstruct() {
    const gauss_rational a ( gauss_rational::integer_type ( 4, 0 ),
                             gauss_rational::integer_type ( 2, 1 ) );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
