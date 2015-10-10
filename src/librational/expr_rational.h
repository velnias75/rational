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
 *
 * @defgroup expr Expression templates
 *
 * The @em expression @em templates @em module allows you to write your code in
 * @em domain-specific @em expressions (DSE) as well as to lazy evaluate them.
 *
 * @b Example: \n
 * To approximate the integral @f$\int_a^b \! \frac{x}{1+x}@f$
 * by evaluating the expression @f$\frac{x}{1+x}@f$
 * for a specified number of equidistant points in the interval @f$\left[ 1, 5\right)@f$,
 * using the @ref gmp, you could write following @c integrate() template function:@code{.cpp}
 * template<class ExprT> inline static
 * Commons::Math::gmp_rational integrate ( const ExprT &e,
 *                                         const Commons::Math::gmp_rational &from,
 *                                         const Commons::Math::gmp_rational &to,
 *                                         std::size_t n ) {
 *     Commons::Math::gmp_rational sum;
 *
 *     const static Commons::Math::gmp_rational two ( 2, 1 );
 *     const Commons::Math::gmp_rational &step ( ( to - from ) / n );
 *
 *     for ( Commons::Math::gmp_rational i ( from + ( step / two ) ); i < to; i += step ) {
 *
 *         sum += Commons::Math::eval_rat_expr ( e, i );
 *     }
 *
 *     return step * sum;
 * }@endcode
 * next you could call the @c integrate() function like: @code{.cpp}
 * // declare the literal 1
 * const Commons::Math::RationalExpressionTraits<Commons::Math::gmp_rational>::expr_type &one (
 *       Commons::Math::mk_rat_lit ( Commons::Math::gmp_rational ( 1, 1 ) ) );
 *
 * // declare the variable x from prototype one
 * const Commons::Math::RationalExpressionTraits<Commons::Math::gmp_rational>::variable_type &x (
 *       Commons::Math::mk_rat_proto_var ( one ) );
 *
 * // approximate the integration for the interval [1, 5) using 10 equidistant points
 * const Commons::Math::gmp_rational &r ( integrate ( x / ( one + x ),
 *                                        Commons::Math::gmp_rational ( 1, 1 ),
 *                                        Commons::Math::gmp_rational ( 5, 1 ), 10 ) );@endcode
 *
 * This would yield the result @f$\frac{422563503196}{145568097675}\approx 2.9@f$ in @c r.
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

    explicit RationalConstant ( const RationalConstant &o ) : c_ ( o.c_ ) {}

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

    RationalBinaryExpression ( const L &l, const H &h );

    explicit RationalBinaryExpression ( const RationalBinaryExpression &o )
        : l_ ( o.l_ ), h_ ( o.h_ ) {}

    ~RationalBinaryExpression();

    result_type operator() ( const Rational<T, GCD, CHKOP> &d );

private:
    typename RationalExpressionTraits<L>::literal_type l_;
    typename RationalExpressionTraits<H>::literal_type h_;
};

template<class T, class L, class H, class OP, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP>
RationalBinaryExpression<T, L, H, OP, GCD, CHKOP>::RationalBinaryExpression ( const L &l,
        const H &h ) : l_ ( l ), h_ ( h ) {}

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

template<class T, class L, class OP, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP>
struct RationalUnaryExpression {

    typedef typename OP::result_type result_type;
    typedef RationalUnaryExpression<T, L, OP, GCD, CHKOP> expr_type;

    RationalUnaryExpression ( const L &l ) : l_ ( l ) {}

    explicit RationalUnaryExpression ( const RationalUnaryExpression &o ) : l_ ( o.l_ ) {}

    ~RationalUnaryExpression() {}

    result_type operator() ( const Rational<T, GCD, CHKOP> &d ) {
        return OP() ( l_ ( d ) );
    }

private:
    typename RationalExpressionTraits<L>::literal_type l_;
};

template<class T, class E, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct RationalExpression {

    typedef typename E::result_type result_type;

    explicit RationalExpression ( const E &e ); // Note: can act as both normal ctor and c-ctor!

    result_type operator() ( const Rational<T, GCD, CHKOP> &d ) {
        return expr_ ( d );
    }

private:
    E expr_;
};

template<class T, class E, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RationalExpression<T, E, GCD, CHKOP>::RationalExpression ( const E &e ) : expr_ ( e ) {}

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

template<class T>
struct _unaryPlus {

    typedef T result_type;

    inline result_type operator() ( const T &d ) const {
        return d;
    }
};

template<class T>
struct _unaryAbs {

    typedef T result_type;

    inline result_type operator() ( const T &d ) const {
        return d.abs();
    }
};

template<class T>
struct _unaryInv {

    typedef T result_type;

    inline result_type operator() ( const T &d ) const {
        return d.inverse();
    }
};

}

}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::Rational<T, GCD, CHKOP>,
        std::plus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
operator+ ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a, const B &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::Rational<T, GCD, CHKOP>,
            std::plus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a,
            Commons::Math::Rational<T, GCD, CHKOP> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::Rational<T, GCD, CHKOP>,
        std::minus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
operator- ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a, const B &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::Rational<T, GCD, CHKOP>,
            std::minus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a,
            Commons::Math::Rational<T, GCD, CHKOP> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::Rational<T, GCD, CHKOP>,
        std::multiplies<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
operator* ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a, const B &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::Rational<T, GCD, CHKOP>,
            std::multiplies<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a,
            Commons::Math::Rational<T, GCD, CHKOP> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
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

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::Rational<T, GCD, CHKOP>,
        std::divides<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
operator/ ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a, const B &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::Rational<T, GCD, CHKOP>,
            std::divides<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a,
            Commons::Math::Rational<T, GCD, CHKOP> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::RationalConstant<B, GCD, CHKOP>,
        std::modulus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator% ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a,
const Commons::Math::Rational<B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::RationalConstant<B, GCD, CHKOP>,
            std::modulus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a,
            Commons::Math::RationalConstant<B, GCD, CHKOP> ( b ) ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalConstant<A, GCD, CHKOP>,
        Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
        std::modulus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator% ( const Commons::Math::Rational<A, GCD, CHKOP> &a,
const Commons::Math::RationalExpression<T, B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalConstant<A, GCD, CHKOP>,
            Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
            std::modulus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT (
                Commons::Math::RationalConstant<A, GCD, CHKOP> ( a ), b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
        std::modulus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
        operator% ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a,
const Commons::Math::RationalExpression<T, B, GCD, CHKOP> &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::RationalExpression<T, B, GCD, CHKOP>,
            std::modulus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a, b ) );
}

template<class T, class A, class B, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalBinaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::Rational<T, GCD, CHKOP>,
        std::modulus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
operator% ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a, const B &b ) {
    typedef Commons::Math::RationalBinaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::Rational<T, GCD, CHKOP>,
            std::modulus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a,
            Commons::Math::Rational<T, GCD, CHKOP> ( b ) ) );
}

template<class T, class A, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalUnaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        std::negate<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
operator- ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a ) {
    typedef Commons::Math::RationalUnaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            std::negate<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a ) );
}

template<class T, class A, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalUnaryExpression<T,
        Commons::Math::RationalConstant<A, GCD, CHKOP>,
        std::negate<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
operator- ( const Commons::Math::Rational<A, GCD, CHKOP> &a ) {
    typedef Commons::Math::RationalUnaryExpression<T,
            Commons::Math::RationalConstant<A, GCD, CHKOP>,
            std::negate<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT (
                Commons::Math::RationalConstant<A, GCD, CHKOP> ( a ) ) );
}

template<class T, class A, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalUnaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::_unaryPlus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
operator+ ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a ) {
    typedef Commons::Math::RationalUnaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::_unaryPlus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a ) );
}

template<class T, class A, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalUnaryExpression<T,
        Commons::Math::RationalConstant<A, GCD, CHKOP>,
        Commons::Math::_unaryPlus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
operator+ ( const Commons::Math::Rational<A, GCD, CHKOP> &a ) {
    typedef Commons::Math::RationalUnaryExpression<T,
            Commons::Math::RationalConstant<A, GCD, CHKOP>,
            Commons::Math::_unaryPlus<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT (
                Commons::Math::RationalConstant<A, GCD, CHKOP> ( a ) ) );
}

template<class T, class A, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalUnaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::_unaryAbs<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
abs ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a ) {
    typedef Commons::Math::RationalUnaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::_unaryAbs<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a ) );
}

template<class T, class A, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalUnaryExpression<T,
        Commons::Math::RationalConstant<A, GCD, CHKOP>,
        Commons::Math::_unaryAbs<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
abs ( const Commons::Math::Rational<A, GCD, CHKOP> &a ) {
    typedef Commons::Math::RationalUnaryExpression<T,
            Commons::Math::RationalConstant<A, GCD, CHKOP>,
            Commons::Math::_unaryAbs<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT (
                Commons::Math::RationalConstant<A, GCD, CHKOP> ( a ) ) );
}

template<class T, class A, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalUnaryExpression<T,
        Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
        Commons::Math::_unaryInv<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
inv ( const Commons::Math::RationalExpression<T, A, GCD, CHKOP> &a ) {
    typedef Commons::Math::RationalUnaryExpression<T,
            Commons::Math::RationalExpression<T, A, GCD, CHKOP>,
            Commons::Math::_unaryInv<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT ( a ) );
}

template<class T, class A, template<typename, bool, template<class, typename, bool> class,
template<typename> class> class GCD, template<class, typename, bool> class CHKOP> inline
Commons::Math::RationalExpression<T, Commons::Math::RationalUnaryExpression<T,
        Commons::Math::RationalConstant<A, GCD, CHKOP>,
        Commons::Math::_unaryInv<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP>, GCD, CHKOP>
inv ( const Commons::Math::Rational<A, GCD, CHKOP> &a ) {
    typedef Commons::Math::RationalUnaryExpression<T,
            Commons::Math::RationalConstant<A, GCD, CHKOP>,
            Commons::Math::_unaryInv<Commons::Math::Rational<T, GCD, CHKOP> >, GCD, CHKOP> ExprT;
    return Commons::Math::RationalExpression<T, ExprT, GCD, CHKOP> ( ExprT (
                Commons::Math::RationalConstant<A, GCD, CHKOP> ( a ) ) );
}

#endif /* COMMONS_MATH_EXPR_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
