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

#define declare_rat_var(type, var) \
    RationalExpression<type, RationalVariable<type> > var ( ( RationalVariable<type>() ) )

namespace Commons {

namespace Math {

template<class ExprT> struct RationalExpressionTraits {
    typedef ExprT expr_type;
    typedef expr_type literal_type;
};

template<class T>
struct RationalVariable {

    typedef Rational<T> result_type;

    result_type operator() ( Rational<T> v ) {
        return v;
    }
};

template<class T>
struct RationalConstant {

    typedef Rational<T> result_type;
    typedef RationalConstant<T> expr_type;

    explicit RationalConstant ( Rational<T> c ) : c_ ( c ) {}

    result_type operator() ( const Rational<T> ) {
        return c_;
    }

private:
    Rational<T> c_;
};

template<class, class> struct RationalExpression;

template<class T, class L, class H, class OP>
struct RationalBinaryExpression {

    typedef typename OP::result_type result_type;
    typedef RationalBinaryExpression<T, L, H, OP> expr_type;

    RationalBinaryExpression ( L l, H h ) : l_ ( l ), h_ ( h ) {}

    ~RationalBinaryExpression();

    result_type operator() ( Rational<T> d );

private:
    typename RationalExpressionTraits<L>::literal_type l_;
    typename RationalExpressionTraits<H>::literal_type h_;
};

template<class T, class L, class H, class OP>
RationalBinaryExpression<T, L, H, OP>::~RationalBinaryExpression() {}

template<class T, class L, class H, class OP>
typename RationalBinaryExpression<T, L, H, OP>::result_type
RationalBinaryExpression<T, L, H, OP>::operator() ( Rational<T> d ) {
    return OP() ( l_ ( d ), h_ ( d ) );
}

template<class T, class E>
struct RationalExpression {

    typedef typename E::result_type result_type;

    explicit RationalExpression ( E e ) : expr_ ( e ) {}

    result_type operator() ( Rational<T> d ) {
        return expr_ ( d );
    }

private:
    E expr_;
};

template<class T> struct RationalExpressionTraits<Rational<T> > {
    typedef RationalConstant<T> literal_type;
    typedef RationalExpression<T, literal_type> expr_type;
};

template<class T>
inline RationalExpression<T, RationalConstant<T> >
mk_rat_lit ( const Rational<T> &r ) {
    return RationalExpression<T, RationalConstant<T> > ( ( RationalConstant<T> ( r ) ) );
}

template<class T> inline
typename RationalExpressionTraits<T>::expr_type::result_type eval_rat_exp ( T expr,
        const typename RationalExpressionTraits<T>::expr_type::result_type &val =
            typename RationalExpressionTraits<T>::expr_type::result_type() ) {
    return ( typename RationalExpressionTraits<T>::expr_type
             ( ( typename RationalExpressionTraits<T>::literal_type ( expr ) ) ) )
           .operator() ( val );
}

}

}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A>, Commons::Math::RationalConstant<B>,
        std::plus<Commons::Math::Rational<T> > > >
operator+ ( const Commons::Math::RationalExpression<T, A> a, Commons::Math::Rational<B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A>, Commons::Math::RationalConstant<B>,
            std::plus<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT ( a,
            Commons::Math::RationalConstant<B> ( b ) ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A>, Commons::Math::RationalExpression<T, B>,
        std::plus<Commons::Math::Rational<T> > > >
operator+ ( Commons::Math::Rational<A> a, const Commons::Math::RationalExpression<T, B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A>, Commons::Math::RationalExpression<T, B>,
            std::plus<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT (
                Commons::Math::RationalConstant<A> ( a ), b ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A>, Commons::Math::RationalExpression<T, B>,
        std::plus<Commons::Math::Rational<T> > > >
operator+ ( Commons::Math::RationalExpression<T, A> a, Commons::Math::RationalExpression<T, B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T, Commons::Math::RationalExpression<T, A>,
            Commons::Math::RationalExpression<T, B>,
            std::plus<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT ( a, b ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A>, Commons::Math::RationalConstant<B>,
        std::minus<Commons::Math::Rational<T> > > >
operator- ( const Commons::Math::RationalExpression<T, A> a, Commons::Math::Rational<B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A>, Commons::Math::RationalConstant<B>,
            std::minus<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT ( a,
            Commons::Math::RationalConstant<B> ( b ) ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A>, Commons::Math::RationalExpression<T, B>,
        std::minus<Commons::Math::Rational<T> > > >
operator- ( Commons::Math::Rational<A> a, const Commons::Math::RationalExpression<T, B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A>, Commons::Math::RationalExpression<T, B>,
            std::minus<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT (
                Commons::Math::RationalConstant<A> ( a ), b ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A>, Commons::Math::RationalExpression<T, B>,
        std::minus<Commons::Math::Rational<T> > > >
operator- ( Commons::Math::RationalExpression<T, A> a, Commons::Math::RationalExpression<T, B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T, Commons::Math::RationalExpression<T, A>,
            Commons::Math::RationalExpression<T, B>,
            std::minus<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT ( a, b ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A>, Commons::Math::RationalConstant<B>,
        std::multiplies<Commons::Math::Rational<T> > > >
operator* ( const Commons::Math::RationalExpression<T, A> a, Commons::Math::Rational<B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A>, Commons::Math::RationalConstant<B>,
            std::multiplies<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT ( a,
            Commons::Math::RationalConstant<B> ( b ) ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A>, Commons::Math::RationalExpression<T, B>,
        std::multiplies<Commons::Math::Rational<T> > > >
operator* ( Commons::Math::Rational<A> a, const Commons::Math::RationalExpression<T, B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A>, Commons::Math::RationalExpression<T, B>,
            std::multiplies<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT (
                Commons::Math::RationalConstant<A> ( a ), b ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A>, Commons::Math::RationalExpression<T, B>,
        std::multiplies<Commons::Math::Rational<T> > > >
operator* ( Commons::Math::RationalExpression<T, A> a, Commons::Math::RationalExpression<T, B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T, Commons::Math::RationalExpression<T, A>,
            Commons::Math::RationalExpression<T, B>,
            std::multiplies<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT ( a, b ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A>, Commons::Math::RationalConstant<B>,
        std::divides<Commons::Math::Rational<T> > > >
operator/ ( const Commons::Math::RationalExpression<T, A> a, Commons::Math::Rational<B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A>, Commons::Math::RationalConstant<B>,
            std::divides<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT ( a,
            Commons::Math::RationalConstant<B> ( b ) ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A>, Commons::Math::RationalExpression<T, B>,
        std::divides<Commons::Math::Rational<T> > > >
operator/ ( Commons::Math::Rational<A> a, const Commons::Math::RationalExpression<T, B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A>, Commons::Math::RationalExpression<T, B>,
            std::divides<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT (
                Commons::Math::RationalConstant<A> ( a ), b ) );
}

template<class T, class A, class B>
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A>, Commons::Math::RationalExpression<T, B>,
        std::divides<Commons::Math::Rational<T> > > >
operator/ ( Commons::Math::RationalExpression<T, A> a, Commons::Math::RationalExpression<T, B> b ) {
    typedef Commons::Math::RationalBinaryExpression<T, Commons::Math::RationalExpression<T, A>,
            Commons::Math::RationalExpression<T, B>,
            std::divides<Commons::Math::Rational<T> > > ExprT;
    return Commons::Math::RationalExpression<T, ExprT> ( ExprT ( a, b ) );
}


#endif /* COMMONS_MATH_EXPR_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
