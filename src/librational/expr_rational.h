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
    typedef expr_type variable_type;
};

template<class T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct RationalVariable {

    typedef Rational<T, GCD, CHKOP> result_type;

    result_type operator() ( const Rational<T, GCD, CHKOP> &v ) const {
        return v;
    }
};

template<class T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct RationalConstant {

    typedef Rational<T, GCD, CHKOP> result_type;

    explicit RationalConstant ( const Rational<T, GCD, CHKOP> &c ) : c_ ( c ) {}

    result_type operator() ( const Rational<T, GCD, CHKOP> & ) const {
        return c_;
    }

private:
    Rational<T, GCD, CHKOP> c_;
};

template<class T, class L, class H, class OP, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP>
struct RationalBinaryExpression {

    typedef typename OP::result_type result_type;
    typedef RationalBinaryExpression<T, L, H, OP, GCD, CHKOP> expr_type;

    RationalBinaryExpression ( const L &l, const H &h ) : l_ ( l ), h_ ( h ) {}

    ~RationalBinaryExpression();

    result_type operator() ( const Rational<T, GCD, CHKOP> &d );

private:
    typename RationalExpressionTraits<L>::literal_type l_;
    typename RationalExpressionTraits<H>::literal_type h_;
};

template<class T, class L, class H, class OP, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP>
RationalBinaryExpression<T, L, H, OP, GCD, CHKOP>::~RationalBinaryExpression() {}

template<class T, class L, class H, class OP, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP>
typename RationalBinaryExpression<T, L, H, OP, GCD, CHKOP>::result_type
RationalBinaryExpression<T, L, H, OP, GCD, CHKOP>::operator() ( const Rational<T, GCD, CHKOP> &d ) {
    return OP() ( l_ ( d ), h_ ( d ) );
}

template<class T, class E, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct RationalExpression {

    typedef typename E::result_type result_type;

    explicit RationalExpression ( const E &e ) : expr_ ( e ) {}

    result_type operator() ( const Rational<T, GCD, CHKOP> &d ) {
        return expr_ ( d );
    }

private:
    E expr_;
};

template<class T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct RationalExpressionTraits<Rational<T, GCD, CHKOP> > {
    typedef RationalConstant<T, GCD, CHKOP> literal_type;
    typedef RationalExpression<T, literal_type, GCD, CHKOP> expr_type;
    typedef RationalExpression<T, RationalVariable<T, GCD, CHKOP>, GCD, CHKOP> variable_type;
};

template<class T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
inline RationalExpression<T, RationalConstant<T, GCD, CHKOP>, GCD, CHKOP>
mk_rat_lit ( const Rational<T, GCD, CHKOP> &r ) {
    return RationalExpression<T, RationalConstant<T, GCD, CHKOP>, GCD, CHKOP>
           ( ( RationalConstant<T, GCD, CHKOP> ( r ) ) );
}

template<class T, class E, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
inline RationalExpression<T, RationalVariable<T, GCD, CHKOP>, GCD, CHKOP>
mk_rat_proto_var ( const RationalExpression<T, E, GCD, CHKOP> & ) {
    return RationalExpression<T, RationalVariable<T, GCD, CHKOP>, GCD, CHKOP>
           ( ( RationalVariable<T, GCD, CHKOP>() ) );
}

template<class T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
inline RationalExpression<T, RationalVariable<T, GCD, CHKOP>, GCD, CHKOP>
mk_rat_proto_var ( const Rational<T, GCD, CHKOP> & ) {
    return RationalExpression<T, RationalVariable<T, GCD, CHKOP>, GCD, CHKOP>
           ( ( RationalVariable<T, GCD, CHKOP>() ) );
}

template<class T> inline
typename RationalExpressionTraits<T>::expr_type::result_type eval_rat_expr ( const T &expr,
        const typename RationalExpressionTraits<T>::expr_type::result_type &val =
            typename RationalExpressionTraits<T>::expr_type::result_type() ) {
    return ( typename RationalExpressionTraits<T>::expr_type
             ( ( typename RationalExpressionTraits<T>::literal_type ( expr ) ) ) )
           .operator() ( val );
}

}

}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::RationalConstant<B, GCD, CHKOP>,
        std::plus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator+ ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a,
const Commons::Math::Rational<B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::RationalConstant<B, GCD, CHKOP>,
            std::plus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a,
            Commons::Math::RationalConstant<B, GCD, CHKOP> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A, GCD, CHKOP>,
        Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
        std::plus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator+ ( const Commons::Math::Rational<A, GCD, CHKOP> &a,
const Commons::Math::RationalExpression<T, B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A, GCD, CHKOP>,
            Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
            std::plus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT (
                Commons::Math::RationalConstant<A, GCD, CHKOP> ( a ), b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
        std::plus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator+ ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a,
const Commons::Math::RationalExpression<T, B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
            std::plus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a, b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::RationalConstant<B, GCD, CHKOP>,
        std::minus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator- ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a,
const Commons::Math::Rational<B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::RationalConstant<B, GCD, CHKOP>,
            std::minus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a,
            Commons::Math::RationalConstant<B, GCD, CHKOP> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A, GCD, CHKOP>,
        Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
        std::minus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator- ( const Commons::Math::Rational<A, GCD, CHKOP> &a,
const Commons::Math::RationalExpression<T, B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A, GCD, CHKOP>,
            Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
            std::minus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT (
                Commons::Math::RationalConstant<A, GCD, CHKOP> ( a ), b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
        std::minus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator- ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a,
const Commons::Math::RationalExpression<T, B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
            std::minus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a, b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::RationalConstant<B, GCD, CHKOP>,
        std::multiplies<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator* ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a,
const Commons::Math::Rational<B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::RationalConstant<B, GCD, CHKOP>,
            std::multiplies<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a,
            Commons::Math::RationalConstant<B, GCD, CHKOP> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A, GCD, CHKOP>,
        Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
        std::multiplies<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator* ( const Commons::Math::Rational<A, GCD, CHKOP> &a,
const Commons::Math::RationalExpression<T, B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A, GCD, CHKOP>,
            Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
            std::multiplies<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT (
                Commons::Math::RationalConstant<A, GCD, CHKOP> ( a ), b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
        std::multiplies<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator* ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a,
const Commons::Math::RationalExpression<T, B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
            std::multiplies<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a, b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::RationalConstant<B, GCD, CHKOP>,
        std::divides<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator/ ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a,
const Commons::Math::Rational<B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::RationalConstant<B, GCD, CHKOP>,
            std::divides<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a,
            Commons::Math::RationalConstant<B, GCD, CHKOP> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A, GCD, CHKOP>,
        Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
        std::divides<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator/ ( const Commons::Math::Rational<A, GCD, CHKOP> &a,
const Commons::Math::RationalExpression<T, B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A, GCD, CHKOP>,
            Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
            std::divides<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT (
                Commons::Math::RationalConstant<A, GCD, CHKOP> ( a ), b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
        std::divides<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator/ ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a,
const Commons::Math::RationalExpression<T, B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
            std::divides<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a, b ) );
}

#endif /* COMMONS_MATH_EXPR_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
