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

/**
 * @file
 * @author Heiko Schäfer <heiko@rangun.de>
 * @copyright 2015 by Heiko Schäfer <heiko@rangun.de>
 */

#ifndef COMMONS_MATH_EXPR_RATIONAL_H
#define COMMONS_MATH_EXPR_RATIONAL_H

#include "rational.h"

namespace Commons {

namespace Math {

template<class ExprT> struct RationalExpressionTraits {
    typedef ExprT expr_type;
    typedef expr_type literal_type;
};

template<class T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD>
struct RationalVariable {

    typedef Rational<T, GCD> result_type;

    result_type operator() ( const Rational<T, GCD> &v ) const {
        return v;
    }
};

template<class T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD>
struct RationalConstant {

    typedef Rational<T, GCD> result_type;
    typedef RationalConstant<T, GCD> expr_type;

    explicit RationalConstant ( const Rational<T, GCD> &c ) : c_ ( c ) {}

    result_type operator() ( const Rational<T, GCD> & ) const {
        return c_;
    }

private:
    Rational<T, GCD> c_;
};

template<class T, class L, class H, class OP, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD>
struct RationalBinaryExpression {

    typedef typename OP::result_type result_type;
    typedef RationalBinaryExpression<T, L, H, OP, GCD> expr_type;

    RationalBinaryExpression ( const L &l, const H &h ) : l_ ( l ), h_ ( h ) {}

    ~RationalBinaryExpression();

    result_type operator() ( const Rational<T, GCD> &d );

private:
    typename RationalExpressionTraits<L>::literal_type l_;
    typename RationalExpressionTraits<H>::literal_type h_;
};

template<class T, class L, class H, class OP, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD>
RationalBinaryExpression<T, L, H, OP, GCD>::~RationalBinaryExpression() {}

template<class T, class L, class H, class OP, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD>
typename RationalBinaryExpression<T, L, H, OP, GCD>::result_type
RationalBinaryExpression<T, L, H, OP, GCD>::operator() ( const Rational<T, GCD> &d ) {
    return OP() ( l_ ( d ), h_ ( d ) );
}

template<class T, class E, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD>
struct RationalExpression {

    typedef typename E::result_type result_type;

    explicit RationalExpression ( const E &e ) : expr_ ( e ) {}

    result_type operator() ( const Rational<T, GCD> &d ) {
        return expr_ ( d );
    }

private:
    E expr_;
};

template<class T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD> struct RationalExpressionTraits<Rational<T, GCD> > {
    typedef RationalConstant<T, GCD> literal_type;
    typedef RationalExpression<T, literal_type, GCD> expr_type;
};

template<class T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD>
inline RationalExpression<T, RationalConstant<T, GCD>, GCD>
mk_rat_lit ( const Rational<T, GCD> &r ) {
    return RationalExpression<T, RationalConstant<T, GCD>, GCD>
           ( ( RationalConstant<T, GCD> ( r ) ) );
}

template<class T, class E, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD>
inline RationalExpression<T, RationalVariable<T, GCD>, GCD>
mk_rat_proto_var ( const RationalExpression<T, E, GCD> & ) {
    return RationalExpression<T, RationalVariable<T, GCD>, GCD> ( ( RationalVariable<T, GCD>() ) );
}

template<class T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD>
inline RationalExpression<T, RationalVariable<T, GCD>, GCD>
mk_rat_proto_var ( const Rational<T, GCD> & ) {
    return RationalExpression<T, RationalVariable<T, GCD>, GCD> ( ( RationalVariable<T, GCD>() ) );
}

template<class T> inline
typename RationalExpressionTraits<T>::expr_type::result_type eval_rat_exp ( const T &expr,
        const typename RationalExpressionTraits<T>::expr_type::result_type &val =
            typename RationalExpressionTraits<T>::expr_type::result_type() ) {
    return ( typename RationalExpressionTraits<T>::expr_type
             ( ( typename RationalExpressionTraits<T>::literal_type ( expr ) ) ) )
           .operator() ( val );
}

}

}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalConstant<B, GCD>,
        std::plus<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator+ ( const Commons::Math::RationalExpression<T, A, GCD> &a,
const Commons::Math::Rational<B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalConstant<B, GCD>,
            std::plus<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT ( a,
            Commons::Math::RationalConstant<B, GCD> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
        std::plus<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator+ ( const Commons::Math::Rational<A, GCD> &a,
const Commons::Math::RationalExpression<T, B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
            std::plus<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT (
                Commons::Math::RationalConstant<A, GCD> ( a ), b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
        std::plus<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator+ ( const Commons::Math::RationalExpression<T, A, GCD> &a,
const Commons::Math::RationalExpression<T, B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T, Commons::Math::RationalExpression<T, A, GCD>,
            Commons::Math::RationalExpression<T, B, GCD>,
            std::plus<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT ( a, b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalConstant<B, GCD>,
        std::minus<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator- ( const Commons::Math::RationalExpression<T, A, GCD> &a,
const Commons::Math::Rational<B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalConstant<B, GCD>,
            std::minus<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT ( a,
            Commons::Math::RationalConstant<B, GCD> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
        std::minus<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator- ( const Commons::Math::Rational<A, GCD> &a,
const Commons::Math::RationalExpression<T, B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
            std::minus<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT (
                Commons::Math::RationalConstant<A, GCD> ( a ), b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
        std::minus<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator- ( const Commons::Math::RationalExpression<T, A, GCD> &a,
const Commons::Math::RationalExpression<T, B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T, Commons::Math::RationalExpression<T, A, GCD>,
            Commons::Math::RationalExpression<T, B, GCD>,
            std::minus<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT ( a, b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalConstant<B, GCD>,
        std::multiplies<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator* ( const Commons::Math::RationalExpression<T, A, GCD> &a,
const Commons::Math::Rational<B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalConstant<B, GCD>,
            std::multiplies<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT ( a,
            Commons::Math::RationalConstant<B, GCD> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
        std::multiplies<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator* ( const Commons::Math::Rational<A, GCD> &a,
const Commons::Math::RationalExpression<T, B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
            std::multiplies<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT (
                Commons::Math::RationalConstant<A, GCD> ( a ), b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
        std::multiplies<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator* ( const Commons::Math::RationalExpression<T, A, GCD> &a,
const Commons::Math::RationalExpression<T, B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T, Commons::Math::RationalExpression<T, A, GCD>,
            Commons::Math::RationalExpression<T, B, GCD>,
            std::multiplies<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT ( a, b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalConstant<B, GCD>,
        std::divides<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator/ ( const Commons::Math::RationalExpression<T, A, GCD> &a,
const Commons::Math::Rational<B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalConstant<B, GCD>,
            std::divides<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT ( a,
            Commons::Math::RationalConstant<B, GCD> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
        std::divides<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator/ ( const Commons::Math::Rational<A, GCD> &a,
const Commons::Math::RationalExpression<T, B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
            std::divides<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT (
                Commons::Math::RationalConstant<A, GCD> ( a ), b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD>, Commons::Math::RationalExpression<T, B, GCD>,
        std::divides<Commons::Math::Rational<T, GCD> >, GCD>, GCD>
        operator/ ( const Commons::Math::RationalExpression<T, A, GCD> &a,
const Commons::Math::RationalExpression<T, B, GCD> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T, Commons::Math::RationalExpression<T, A, GCD>,
            Commons::Math::RationalExpression<T, B, GCD>,
            std::divides<Commons::Math::Rational<T, GCD> >, GCD> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD> ( ExprT ( a, b ) );
}


#endif /* COMMONS_MATH_EXPR_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
