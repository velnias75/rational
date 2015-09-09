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

#include <cppunit/ui/text/TestRunner.h>

#include "rationaltest.h"

int main ( int, char ** ) {

    CppUnit::TextUi::TestRunner runner;
    CppUnit::TestFactoryRegistry &registry = CppUnit::TestFactoryRegistry::getRegistry();

    runner.addTest ( registry.makeTest() );

    return !runner.run ( "", false );
}

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
