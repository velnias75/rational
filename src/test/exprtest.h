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

#include <cppunit/extensions/HelperMacros.h>

#include "gmp_rational.h"

namespace Commons {

namespace Math {

template<class T>
struct RationalVar {

    Rational<T> &operator () ( Rational<T> &v ) {
        return v;
    }
};

template<class T>
struct RationalConstant {

    explicit RationalConstant ( const Rational<T> &r ) : c ( r ) {}

    Rational<T> operator () ( const Rational<T> & ) {
        return c;
    }

private:
    Rational<T> c;
};

template<class T, class L, class H, class OP>
struct RationalBinaryExpression {

    RationalBinaryExpression ( const L &l, const H &h ) : l_ ( l ), h_ ( h ) {}

    ~RationalBinaryExpression();

    Rational<T> operator () ( const Rational<T> &d );

private:
    L l_;
    H h_;
};

template<class T, class L, class H, class OP>
RationalBinaryExpression<T, L, H, OP>::~RationalBinaryExpression() {}

template<class T, class L, class H, class OP>
Rational<T> RationalBinaryExpression<T, L, H, OP>::operator () ( const Rational<T> &d ) {
    return OP() ( l_ ( d ), h_ ( d ) );
}

template<class T, class E>
struct RationalExpression {

    explicit RationalExpression ( const E &e ) : expr_ ( e ) {}

    Rational<T> operator() ( const Rational<T> &d = Rational<T>() ) {
        return expr_ ( d );
    }

private:
    E expr_;
};

template<class T>
inline RationalExpression<T, RationalConstant<T> > mk_rat_lit ( const Rational<T> &r ) {
    return RationalExpression<T, RationalConstant<T> > ( ( RationalConstant<T> ( r ) ) );
}

template<class T, class N>
inline RationalExpression<T, RationalConstant<T> > mk_rat_lit ( const N &n, const N &d ) {
    return RationalExpression<T, RationalConstant<T> > ( ( RationalConstant<T>
            ( Rational<T> ( n, d ) ) ) );
}

template<class T, class N> inline
RationalExpression<T, RationalConstant<T> > mk_rat_lit ( const N &w, const N &n, const N &d ) {
    return RationalExpression<T, RationalConstant<T> > ( ( RationalConstant<T>
            ( Rational<T> ( w, n, d ) ) ) );
}

}

}

template<class T, class A, class B>
Commons::Math::RationalExpression<T,
        Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A>,
        Commons::Math::RationalExpression<T, B>,
        std::plus<Commons::Math::Rational<T> > > >
        operator+ ( const Commons::Math::RationalExpression<T, A> &a,
const Commons::Math::RationalExpression<T, B> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T, Commons::Math::RationalExpression<T, A>,
            Commons::Math::RationalExpression<T, B>, std::plus<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT ( a, b ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T,
        Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A>,
        Commons::Math::RationalExpression<T, B>,
        std::minus<Commons::Math::Rational<T> > > >
        operator- ( const Commons::Math::RationalExpression<T, A> &a,
const Commons::Math::RationalExpression<T, B> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T, Commons::Math::RationalExpression<T, A>,
            Commons::Math::RationalExpression<T, B>, std::minus<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT ( a, b ) );
}

#pragma GCC diagnostic ignored "-Winline"
#pragma GCC diagnostic ignored "-Weffc++"
#pragma GCC diagnostic push
class ExprTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE ( ExprTest );
    CPPUNIT_TEST ( testExpression );
    CPPUNIT_TEST_SUITE_END();

public:
    ExprTest();

    void setUp();
    void tearDown();

    void testExpression();
};
#pragma GCC diagnostic pop

#endif /* EXPRTESTCASE_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
