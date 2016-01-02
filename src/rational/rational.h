/*
 * Copyright 2015-2016 by Heiko Schäfer <heiko@rangun.de>
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
 * @copyright 2015-2016 by Heiko Schäfer <heiko@rangun.de>
 *
 * @defgroup main Rational fraction class
 *
 * Include `rational.h` to be able to do fraction calculations. By simply including `rational.h`
 * and specifying the storage type (any integer variant) you can create and use a fractional data
 * type. For example, `Commons::Math::Rational<long> foo(3, 4)` would create a fraction named `foo`
 * with a value of @f$ \frac{3}{4} @f$, and store the fraction using the `long` data type.\n \n
 *
 * The @em storage @em type should represent all integers within some (possibly infinite) interval
 * @f$[0, 1]@f$. For example, the native @c signed or @c unsigned @c int and @c long types, or
 * arbitrary-precision integers, may be used. Beyond ordinary integers, you should also be
 * able to use any other [Euclidean domain](https://en.wikipedia.org/wiki/Euclidean_domain),
 * perhaps not even an [ordered ring](https://en.wikipedia.org/wiki/Ordered_ring), but support
 * for such types is experimental and has not been thoroughly tested. In fact, you should be able
 * to use any [integral domain](https://en.wikipedia.org/wiki/Integral_domain), but you may need
 * to apply a more sophisticated GCD algorithm (see also @ref gcd); you can fall back to
 * Commons::Math::GCD_null if overflow is not a concern in practice. Finally, using non-integral
 * domains is very likely to fail.\n \n
 *
 * @b Example: \n To approximate the @b reciprocal of the @em golden @em ratio
 * (@f$ \varphi = \phi^{-1} @f$) \n @f$ \varphi = \frac{\sqrt{5}-1}{2} @f$ \n by iteratively
 * calculating \n @f$x_{n}=\frac{F_{n+1}}{F_{n}}@f$, where @f$F_{n}@f$ is the @em n-th
 * Fibonacci number\n you could write: @code{.cpp}
 * Rational<uint64_t> phi ( 1u, 1u ); // init with any Fibonacci(n+1), Fibonacci(n)
 *
 * for ( std::size_t i = 0u; i < 91u; ++i ) ( ++phi ).invert(); // i >= 91u exceeds uint64@endcode
 * which will result in @f$ \varphi \approx
 * \frac{7540113804746346429}{12200160415121876738} = 0.61803398874989484820458683436563811772 @f$
 *
 * @note Use @c Commons::Math::Rational::invert() or just add @c 1 to get @f$ \phi @f$
 */

#ifndef COMMONS_MATH_RATIONAL_H
#define COMMONS_MATH_RATIONAL_H

#include <functional>
#include <algorithm>
#include <sstream>
#include <cstdlib>
#include <string>
#include <limits>
#include <vector>
#include <stack>
#include <cmath>

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
#include <type_traits>
#endif

#ifdef __EXCEPTIONS
#include <stdexcept>
#endif

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
#define RATIONAL_NOEXCEPT  noexcept
#define RATIONAL_CONSTEXPR constexpr
#else
#define RATIONAL_NOEXCEPT  throw()
#define RATIONAL_CONSTEXPR
#endif

namespace Commons {

namespace tmp {

#pragma GCC diagnostic ignored "-Wctor-dtor-privacy"
#pragma GCC diagnostic push
template<typename T>
class _isClassT {
private:
    typedef char One;
    typedef struct {
        char a[2];
    } Two;
    template<typename C> static One test ( int C:: * );
    template<typename C> static Two test ( ... );
public:
    enum { Yes = sizeof ( _isClassT<T>::template test<T> ( 0L ) ) == 1 };
    enum { No = !Yes };
};
#pragma GCC diagnostic pop

template<bool C, typename Ta, typename Tb>
struct _ifThenElse;

template<typename Ta, typename Tb>
struct _ifThenElse<true, Ta, Tb> {
    typedef Ta ResultT;
};

template<typename Ta, typename Tb>
struct _ifThenElse<false, Ta, Tb> {
    typedef Tb ResultT;
};

}

namespace Math {

template<typename, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class, template<class, typename, bool> class, bool> struct _pow;

template<typename, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class, template<class, typename, bool> class, bool> struct _mod;

template<typename, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class, template<class, typename, bool> class, bool> struct _abs;

template<typename, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class, template<class, typename, bool> class, bool> struct _lcm;

template<typename, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class, template<class, typename, bool> class, bool>
struct _swapSign;

/**
 * @ingroup main
 * @brief Type coversion policy class
 *
 * Specialize this class to convert a type from @c T to @c U
 *
 * @tparam T the type to convert from
 */
template<typename T> struct TYPE_CONVERT {

    /**
     * @brief Constructs a type converter
     *
     * @param[in] v the value to convert
     */
    inline explicit TYPE_CONVERT ( const T& v ) : val ( v ) {}

    /**
     * @brief converts the value to @c U
     *
     * @tparam U the type to convert to
     */
    template<typename U> inline U convert() const {
        return static_cast<U> ( val );
    }

private:
    const T& val;
};

template<> struct TYPE_CONVERT<std::string> {

    /**
     * @brief Constructs a type converter
     *
     * @param[in] v the value to convert
     */
    inline explicit TYPE_CONVERT ( const std::string &v ) : val ( v ) {}

    /**
     * @brief converts the value to @c U
     *
     * @tparam U the type to convert to
     */
    template<typename U> inline U convert() const {
        U aux;
        ( std::istringstream ( val ) ) >> aux;
        return aux;
    }

private:
    const std::string &val;
};

template<> struct TYPE_CONVERT<const char *> {

    /**
     * @brief Constructs a type converter
     *
     * @param[in] v the value to convert
     */
    inline explicit TYPE_CONVERT ( const char *f, const char *l = 0L ) :
        val ( f ), len ( l ) {}

    /**
     * @brief converts the value to @c U
     *
     * @tparam U the type to convert to
     */
    template<typename U> inline U convert() const {
        U aux;
        ( std::istringstream ( len ? std::string ( val, len ) : val ) ) >> aux;
        return aux;
    }

private:
    const char *val;
    const char *len;
};

template<> struct TYPE_CONVERT<float> {

    inline explicit TYPE_CONVERT ( const float v ) : val ( v ) {}

    template<class U> inline U convert() const {
        return static_cast<U> ( val );
    }

private:
    const float val;
};

template<> struct TYPE_CONVERT<double> {

    inline explicit TYPE_CONVERT ( const double v ) : val ( v ) {}

    template<class U> inline U convert() const {
        return static_cast<U> ( val );
    }

private:
    const double val;
};

template<> struct TYPE_CONVERT<long double> {

    inline explicit TYPE_CONVERT ( const long double v ) : val ( v ) {}

    template<class U> inline U convert() const {
        return static_cast<U> ( val );
    }

private:
    const long double val;
};

/**
 * @ingroup main
 * @brief @c %EPSILON for float approximation
 *
 * Specialize this class if you need another @c %EPSILON (error tolerance)
 *
 * By default this is @code std::numeric_limits<T>::epsilon() @endcode
 *
 * @tparam T storage type of Rational
 */
template<typename T> struct EPSILON {

    /**
     * @brief the value of @c %EPSILON
     */
    inline static const T value() {
        return std::numeric_limits<T>::epsilon();
    }
};

template<typename R> struct SQRT_HERON_ITERATE {

    inline bool operator() ( const typename R::integer_type &x,
                             const typename R::integer_type &y ) const {
        return ! ( x > std::numeric_limits<typename R::integer_type>::max() / y );
    }

    inline bool operator() ( const R &x, const R &y ) const {

        R v ( x );

        return ! ( ( v.numerator() > std::numeric_limits<typename R::integer_type>::max() /
                     y.denominator() ) ||
                   ( v.denominator() > std::numeric_limits<typename R::integer_type>::max() /
                     y.numerator() ) ||
                   ( std::numeric_limits<typename R::integer_type>::max() -
                     ( v.numerator() * y.denominator() ) < ( v.denominator() * y.numerator() ) ) ||
                   ( ( v += y ).denominator() >
                     std::numeric_limits<typename R::integer_type>::max() / two_ ) );
    }

private:
    static typename R::integer_type two_;
};

template<typename R> typename R::integer_type SQRT_HERON_ITERATE<R>::two_ ( R::one_ + R::one_ );

template<typename, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class, template<class, typename, bool> class, typename, bool,
         template<typename> class = EPSILON,
         template<typename> class = TYPE_CONVERT> struct _approxFract;
/**
 * @ingroup main
 * @brief unchecked operator
 *
 * Delegates the operator @c Op without any overflow/wrap check
 *
 * @tparam Op operator functor
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 */
template<class Op, typename T, bool IsSigned>
struct NO_OPERATOR_CHECK {

    RATIONAL_CONSTEXPR inline T operator() ( const T &x, const T& y ) const {
        return Op() ( x, y );
    }
};

template<typename T, bool IsSigned>
struct NO_OPERATOR_CHECK<std::negate<T>, T, IsSigned> {

    RATIONAL_CONSTEXPR inline T operator() ( const T &x ) const {
        return std::negate<T>() ( x );
    }
};

/**
 * @ingroup main
 * @brief checked operator
 *
 * Checks the operands on signed overflows, resp. unsigned wraps
 * and throws a @c std::domain_error if the check fails, else
 * delegates to the operator @c Op
 *
 * @tparam Op operator functor
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 */
template<class Op, typename T, bool IsSigned>
struct ENABLE_OVERFLOW_CHECK {

    RATIONAL_CONSTEXPR inline T operator() ( const T &x, const T& y ) const {
        return NO_OPERATOR_CHECK<Op, T, IsSigned>() ( x, y );
    }
};

/**
 * @ingroup main
 * @ingroup gcd
 * @brief NULL GCD algorithm implementation
 *
 * Despite it's name this GCD does uncondionally return @c T(1)
 *
 * It can be useful if reduction of fractions is not wanted.
 *
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, bool IsSigned, template<class, typename = T, bool = IsSigned> class CHKOP,
         template<typename> class CONV = TYPE_CONVERT> struct GCD_null;

/**
 * @ingroup main
 * @ingroup gcd
 * @brief Stein GCD algorithm implementation
 *
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, bool IsSigned, template<class, typename = T, bool = IsSigned> class CHKOP,
         template<typename> class CONV = TYPE_CONVERT> struct GCD_stein;

/**
 * @ingroup main
 * @ingroup gcd
 * @brief Euclid GCD algorithm (safe) implementation
 *
 * This implementation supports overlow/wrap checking
 *
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, bool IsSigned, template<class, typename = T, bool = IsSigned> class CHKOP,
         template<typename> class CONV = TYPE_CONVERT> struct GCD_euclid;

/**
 * @ingroup main
 * @ingroup gcd
 * @brief Euclid GCD algorithm (fast) implementation
 *
 * @see GCD_euclid if your number class doesn't support all needed operators
 *
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, bool IsSigned, template<class, typename = T, bool = IsSigned> class CHKOP,
         template<typename> class CONV = TYPE_CONVERT> struct GCD_euclid_fast;

/**
 * @ingroup main
 * @brief Traits struct to choose NumberType for expression evaluation
 *
 * If this trait isn't specialized, it defaults to @c long @c double
 *
 * @tparam T integer_type to get corresponding @c NumberType for
 */
template<typename T>
struct ExpressionEvalTraits {
    typedef long double NumberType; ///< the corresponding @c NumberType
};

template<typename T> struct _type_round_helper {
    inline typename ExpressionEvalTraits<T>::NumberType
    operator() ( const typename ExpressionEvalTraits<T>::NumberType &tr ) const {
        return typename ExpressionEvalTraits<T>::NumberType ( 0.5 ) + tr;
    }
};

/**
 * @ingroup main
 * @brief %Rational (fraction) template class
 *
 * @note All %Rational objects are reduced (see @ref gcd)
 *
 * @tparam T storage type
 * @tparam GCD GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD = GCD_euclid_fast,
         template<class, typename = T, bool = std::numeric_limits<T>::is_signed>
         class CHKOP = NO_OPERATOR_CHECK>
class Rational {

#if defined(__GXX_EXPERIMENTAL_CXX0X__) || __cplusplus >= 201103L
    static_assert ( std::is_integral<T>::value || std::numeric_limits<T>::is_integer,
                    "only integer types are allowed as storage type" );
#endif

    friend struct _swapSign<T, GCD, CHKOP, std::numeric_limits<T>::is_signed>;
    friend struct _mod<T, GCD, CHKOP, std::numeric_limits<T>::is_signed>;
    template<typename, template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class, template<class, typename, bool> class, typename, bool,
             template<typename> class, template<typename> class> friend struct _approxFract;
public:
    /**
     * @brief storage type
     */
    typedef T integer_type;

    /**
     * @brief type of the return value of mod()
     *
     * This type is based on @c std::pair, where @c first is the integral value and @c second the
     * fractional part
     */
    typedef typename _mod<integer_type, GCD, CHKOP,
            std::numeric_limits<integer_type>::is_signed>::pair_type mod_type;

    /**
     * @brief addition (@c a @b + @c b) operator
     */
    typedef CHKOP<std::plus<T> > op_plus;

    /**
     * @brief subtraction (@c a @b - @c b) operator
     */
    typedef CHKOP<std::minus<T> > op_minus;

    /**
     * @brief subtraction (@c a @b - @c b) operator
     */
    typedef CHKOP<std::negate<T> > op_negate;

    /**
     * @brief multiplication (@c a @b * @c b) operator
     */
    typedef CHKOP<std::multiplies<T> > op_multiplies;

    /**
     * @brief division (@c a @b / @c b) operator
     */
    typedef CHKOP<std::divides<T> > op_divides;

    /**
     * @brief modulus (@c a @b % @c b) operator
     */
    typedef CHKOP<std::modulus<T> > op_modulus;

    /**
     * @brief creates a default (null) %Rational
     *
     * Creates a fraction @f$ \frac{0}{1} @f$
     */
    RATIONAL_CONSTEXPR Rational() : m_numer (), m_denom ( one_ ) {}

    /**
     * @brief creates a %Rational
     *
     * Creates a copy of @c numer and divides it by @c denom
     *
     * @note to achieve a continued faction like \f$\frac{1}{\frac{1}{2}}\f$ you must explicitely
     * cast the @em numerator, i.e.
     * @code const Rational<rational_type> x ( Rational<rational_type> ( 1 ),
     *   Rational<rational_type> ( 1,2 )); @endcode
     *
     * @tparam U1 GCD algorithm of the numerator
     * @tparam V1 operator checker of the numerator
     *
     * @tparam U2 GCD algorithm of the denominator
     * @tparam V2 operator checker of the denominator
     *
     * @param[in] numer the numerator
     * @param[in] denom the denominator
     */
    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U1,
             template<class, typename, bool> class V1,
             template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U2,
             template<class, typename, bool> class V2>
    Rational ( const Rational<integer_type, U1, V1> &numer,
               const Rational<integer_type, U2, V2> &denom ) :
        m_numer ( numer.numerator() ), m_denom ( numer.denominator() ) {
        *this *= denom.inverse();
    }

    /**
     * @brief creates a %Rational
     *
     * Creates a fraction @f$ \frac{numer}{denom} @f$
     *
     * @param[in] numer the numerator
     * @param[in] denom the denominator
     */
    Rational ( const integer_type &numer, const integer_type &denom );

    /**
     * @brief creates a inproper (mixed) %Rational
     *
     * Creates a fraction @f$ whole + \frac{numer}{denom} @f$
     *
     * @param[in] whole whole number part
     * @param[in] numer the numerator
     * @param[in] denom the denominator
     */
    Rational ( const integer_type &whole, const integer_type &numer,
               const integer_type &denom ) : m_numer (), m_denom ( one_ ) {
        *this += Rational ( numer, denom ) += Rational ( whole );
    }

    /**
     * @brief creates an approximated %Rational
     *
     * @b Example: \n @c Commons::Math::Rational<uint64_t>( std::sqrt( 2.0l ) )
     * to get @f$ \sqrt{2} \approx \frac{6333631924}{4478554083} @f$
     *
     * @see Commons::Math::EPSILON to control the quality of the approximation
     *
     * @tparam NumberType type of the number to approximate
     *
     * @param[in] number the number to create an approximated %Rational of
     */
    template<typename NumberType>
    Rational ( const NumberType &number );

    /**
     * @brief creates a %Rational approximated by an expression
     *
     * If @c expr is not @c null and not empty then it is parsed and evaluated as
     * @c long @c double expression and approximated to a fraction.
     *
     * The @em expression can be a simple infix arithmetic expression containing
     *
     * * addition (@c +); also @em unary
     * * subtraction (@c -); also @em unary
     * * multiplication (@c *)
     * * division (@c /)
     * * modulus (@c %)
     * * parenthesises
     *
     * Numbers can be integers or floats in non-scientific notation. Allowed are spaces, tabs
     * and newlines around numbers, parenthesises and operators.
     *
     * The expression gets evaluated into a sequence of Commons::Math::Rational terms. Floats
     * are approximated using @c ExpressionEvalTraits<integer_type> as float number type.\n
     * @c ExpressionEvalTraits<integer_type> corresponds to @c long @c double if not specialized
     *
     * In case of errors an @c std::runtime_exception is thrown if exceptions are enabled,
     * else the result is undefined.
     *
     * @b Example: \n @code{.cpp}
     * Rational<long> x ( "(11/2) * +(4.25+3.75)" );@endcode produces the fraction
     * @f$x = \frac{44}{1}@f$
     *
     * @see Rational(const NumberType &number)
     * @see ExpressionEvalTraits
     *
     * @param[in] expr the expression to evaluate and approximate
     */
    Rational ( const char *expr );

    ~Rational();

    /**
     * @brief assigns from a @c NumberType
     *
     * The @c number is approximated to a %Rational, then it gets assigned
     *
     * @see Rational(const NumberType &number)
     * @see Commons::Math::EPSILON
     *
     * @tparam NumberType type of the number to approximate
     *
     * @param[in] number the number to assign
     */
    template<typename NumberType>
    inline Rational &operator= ( const NumberType &number ) {
        return ( *this = Rational ( number ) );
    }

#ifndef __clang__
    inline operator void() const RATIONAL_NOEXCEPT {}
#endif

    /**
     * @brief convert to @c NumberType
     *
     * @tparam NumberType type of the number to approximate
     *
     * @return the number value of the %Rational
     */
    template<typename NumberType>
    RATIONAL_CONSTEXPR inline operator NumberType() const {
        return TYPE_CONVERT<integer_type> ( m_numer ).template convert<NumberType>() /
               TYPE_CONVERT<integer_type> ( m_denom ).template convert<NumberType>();
    }

    /**
     * @brief gets the numerator
     *
     * @return the numerator
     */
    RATIONAL_CONSTEXPR inline integer_type numerator() const RATIONAL_NOEXCEPT {
        return m_numer;
    }

    /**
     * @brief gets the denominator
     *
     * @return the denominator
     */
    RATIONAL_CONSTEXPR inline integer_type denominator() const RATIONAL_NOEXCEPT {
        return m_denom;
    }

    /**
     * @brief Structure holding a description of a repeating fraction
     *
     * @see Commons::Math::Rational::decompose()
     * @see Commons::Math::Rational::rf()
     */
    typedef struct _rf_info {

        inline _rf_info ( const integer_type &r, std::size_t lz = 0u,
                          const integer_type &p = integer_type(),  std::size_t plz = 0u ) :
            reptend ( r ), leading_zeros ( lz ), pre ( p ), pre_leading_zeros ( plz ),
            pre_digits(), reptend_digits() {}

        inline _rf_info() : reptend(), leading_zeros ( 0u ), pre(), pre_leading_zeros ( 0u ),
            pre_digits(), reptend_digits() {}

        integer_type reptend; ///< the repeating part as integer_type
        std::size_t leading_zeros; ///< the amount of zeros at the beginning of @c reptend
        integer_type pre; ///< the digits before @c reptend as integer_type
        std::size_t pre_leading_zeros; ///< the amount of zeros at the beginning of @c pre

        std::vector<integer_type> pre_digits; ///< the part before the reptent as digit sequence
        std::vector<integer_type> reptend_digits; ///< the repeating part as digit sequence

    } rf_info;

    /**
     * @ingroup main
     * @brief Constructs a fraction from a repeating decimal
     *
     * @see Commons::Math::Rational::decompose()
     *
     * The fraction is calculated by the formula: \n \n
     * @f$ \frac{\displaystyle{\mathrm{pre}} + \frac{\displaystyle{\mathrm{reptend}}}{\begin{cases}
     * \displaystyle{1} & \displaystyle{\text{if } \mathrm{x} = 0} \\
     * \displaystyle{10^{\displaystyle{\lceil\log_{10}(|\mathrm{reptend}| + 1)\rceil +
     * \mathrm{leading\_zeros}}} - 1} & \displaystyle{\text{if } \mathrm{x} \neq 0}
     * \end{cases}}}{\displaystyle{10}^{\displaystyle{\displaystyle{\lceil\log_{10}(|\mathrm{pre}|
     * + 1)\rceil + \mathrm{pre\_leading\_zeros}}}}} @f$
     *
     * @remarks
     * * to get an intuitive result @c reptend and @c pre should be positive numbers
     * * the resulting fraction will be within @f$ 0 \leq x \leq 1@f$, where @f$ x @f$ is the
     * decimal value of the fraction
     *
     * @b Examples: \n
     * * to construct a fraction representing
     *   @f$\frac{13717421}{111111111} = 0.\overline{123456789}@f$ you'll need to write: @code{.cpp}
     * Commons::Math::Rational<long> frac =
     *      Commons::Math::Rational<long>(Commons::Math::rf_info(123456789));@endcode
     * * to construct a fraction representing
     *   @f$\frac{667}{6000} = 0.1111\overline{66}@f$ you'll need to write: @code{.cpp}
     * Commons::Math::Rational<long> frac =
     *      Commons::Math::Rational<long>(Commons::Math::rf_info(6, 0, 1111));@endcode
     *
     * @param[in] info a repeating decimal description
     */
    Rational ( const rf_info &info );

    /**
     * @brief Splits a fraction in its whole and repetitive part
     *
     * @param[out] rf_info rf_info structure to store the result
     * @param[in] base the base, defaults to integer_type (10) as a decimal base
     *
     * @return the whole part of the fraction
     */
    integer_type decompose ( rf_info &rf_info,
                             const integer_type &base = integer_type ( 10 ) ) const;
    /**
     * @brief extract the integral and fractional part
     *
     * each part has the same sign as the Rational
     *
     * @see mod_type
     *
     * @return a @c mod_type containing the integral and fractional part
     */
    RATIONAL_CONSTEXPR mod_type mod() const;

    /**
     * @brief gets the absolute %Rational
     *
     * Returns for @em signed types
     * @c ( this->numerator() < T() ? -(*this) : (*this) )
     * and for @em unsigned types @c (*this)
     *
     * @return a copy of the absolute %Rational
     */
    RATIONAL_CONSTEXPR inline Rational abs() const {
        return _abs<integer_type, GCD, CHKOP,
               std::numeric_limits<integer_type>::is_signed>() ( *this );
    }

    /**
     * @brief gets the a new %Rational raised to the power of @c exp
     *
     * If overflows can get ruled out (i.e. with use of Commons::Math::gmp_rational),
     * you can speed up the calculation by initializing a %Rational with
     * Commons::Math::GCD_null and initializing the resulting power in a %Rational
     * with the desired GCD: @code{.cpp}
     * const Commons::Math::gmp_rational r(3, 4);
     *
     * typedef Commons::Math::Rational<Commons::Math::gmp_rational::integer_type,
     *         Commons::Math::GCD_null, Commons::Math::NO_OPERATOR_CHECK> gmp_nogcd_rational;
     *
     * const gmp_nogcd_rational &nr(gmp_nogcd_rational(r.numerator(), r.denominator()).pow(256));
     *
     * std::cout << Commons::Math::gmp_rational(nr.numerator(), nr.denominator())
     *           << std::endl;@endcode
     *
     * @warning negative values or zero can cause undefined behaviour if compiled
     * without exceptions
     *
     * @param[in] exp the exponent to raise this %Rational
     * @return a copy of the absolute %Rational raised to the power of @c exp
     */
    inline Rational pow ( const integer_type &exp ) const {
        return _pow<integer_type, GCD, CHKOP,
               std::numeric_limits<integer_type>::is_signed>() ( *this, exp );
    }

    /**
     * @brief gets the square root as a new %Rational
     *
     * @return a square root copy of %Rational
     */
    inline Rational sqrt() const {

        if ( m_numer == m_denom ) return *this;

        const Rational half ( one_, one_ + one_ );
        Rational aux, inv, x ( ( Rational ( one_, one_ ) += *this ) *= half );

        while ( SQRT_HERON_ITERATE<Rational>() ( m_numer, ( inv = x.inverse() ).denominator() ) &&
                SQRT_HERON_ITERATE<Rational>() ( x, ( aux = Rational ( *this ) *= inv ) ) ) {

            x += aux;
            x *= half;
        }

        typename tmp::_ifThenElse<tmp::_isClassT<typename mod_type::first_type>::Yes,
                 const typename mod_type::first_type &,
                 const typename mod_type::first_type>::ResultT m ( x.mod().first );

        if ( m != typename mod_type::first_type() ) {

            const Rational &y ( m );

            return ( y * y ) == *this ? y : x;
        }

        return x;
    }

    /**
     * @brief inverts the %Rational
     *
     * @b Example: \n @f$ \frac{-5}{12} @f$ will become to @f$ \frac{-12}{5} @f$
     *
     * @return the inverted %Rational
     */
    Rational &invert();

    /**
     * @brief gets a copy of the inverted %Rational
     *
     * @see invert()
     *
     * @return a copy of the inverted %Rational
     */
    RATIONAL_CONSTEXPR Rational inverse() const;

    /**
     * @brief add and assign a %Rational
     *
     * @param[in] other the %Rational to add and assign
     *
     * @return the %Rational
     */
    inline Rational& operator+= ( const Rational& other ) {
        return knuth_addSub<op_plus> ( other );
    }

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend Rational &operator+= ( Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( x += Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief add a %Rational
     *
     * @param[in] other the %Rational to add
     *
     * @return a new %Rational
     */
    RATIONAL_CONSTEXPR Rational operator+ ( const Rational& other ) const;

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend Rational operator+ ( const Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( Rational ( x ) += Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief get a copy of the %Rational
     *
     * @return a copy of %Rational
     */
    RATIONAL_CONSTEXPR inline Rational operator+() const {
        return Rational ( *this );
    }

    /**
     * @brief pre-increment a %Rational
     *
     * the result will be
     * @f$ \frac{numerator + denominator}{denominator} \Rightarrow numerator + 1@f$
     *
     * @return the incremented %Rational
     */
    inline Rational& operator++() {
        m_numer = op_plus() ( m_numer, m_denom );
        return reduce();
    }

    /**
     * @brief post-increment a %Rational
     *
     * @see operator++()
     *
     * @return a copy of %Rational
     */
    inline Rational operator++ ( int ) {
        Rational tmp ( *this );
        ++*this;
        return tmp;
    }

    /**
     * @brief subtract and assign a %Rational
     *
     * @param[in] other the %Rational to subtract and assign
     *
     * @return the %Rational
     */
    inline Rational& operator-= ( const Rational& other ) {
        return knuth_addSub<op_minus> ( other );
    }

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend Rational &operator-= ( Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( x -= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief subtract a %Rational
     *
     * @param[in] other the %Rational to subtract
     *
     * @return a new %Rational
     */
    RATIONAL_CONSTEXPR inline Rational operator- ( const Rational& other ) const {
        return ( Rational ( *this ) -= other );
    }

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend Rational operator- ( const Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( Rational ( x ) -= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief get a negated copy of the %Rational
     *
     * @return a negated copy of %Rational
     */
    inline Rational operator-() const {
        Rational tmp ( *this );
        tmp.m_numer = op_negate() ( tmp.m_numer );
        return tmp;
    }

    /**
      * @brief pre-decrement a %Rational
      *
      * the result will be
      * @f$ \frac{numerator - denominator}{denominator} \Rightarrow numerator - 1@f$
      *
      * @return the incremented %Rational
      */
    inline Rational& operator--() {
        m_numer = op_minus() ( m_numer, m_denom );
        return reduce();
    }

    /**
     * @brief post-decrement a %Rational
     *
     * @see operator--()
     *
     * @return a copy of %Rational
     */
    inline Rational operator-- ( int ) {
        Rational tmp ( *this );
        --*this;
        return tmp;
    }

    /**
     * @brief multiply and assign a %Rational
     *
     * @param[in] other the %Rational to multiply and assign
     *
     * @return the %Rational
     */
    Rational& operator*= ( const Rational& other );

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend Rational &operator*= ( Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( x *= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief multiply a %Rational
     *
     * @param[in] other the %Rational to multiply
     *
     * @return a new %Rational
     */
    RATIONAL_CONSTEXPR Rational operator* ( const Rational& other ) const;

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend Rational operator* ( const Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( Rational ( x ) *= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief divide and assign a %Rational
     *
     * @param[in] other the %Rational to divide and assign
     *
     * @return the %Rational
     */
    inline Rational& operator/= ( const Rational& other ) {
        return ( *this *= other.inverse() );
    }

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend Rational &operator/= ( Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( x /= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief divide a %Rational
     *
     * @param[in] other the %Rational to divide
     *
     * @return a new %Rational
     */
    RATIONAL_CONSTEXPR inline Rational operator/ ( const Rational& other ) const {
        return ( Rational ( *this ) /= other );
    }

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend Rational operator/ ( const Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( Rational ( x ) /= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief modulo and assign a %Rational
     *
     * @param[in] other the %Rational to modulo and assign
     *
     * @return the %Rational
     */
    Rational& operator%= ( const Rational& other );

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend Rational &operator%= ( Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( x %= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief modulo a %Rational
     *
     * @param[in] other the %Rational to modulo
     *
     * @return a new %Rational
     */
    RATIONAL_CONSTEXPR inline Rational operator% ( const Rational& other ) const {
        return ( Rational ( *this ) %= other );
    }

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend Rational operator% ( const Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( Rational ( x ) %= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test on equality
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if equal, @c false otherwise
     */
    RATIONAL_CONSTEXPR inline bool operator== ( const Rational &other ) const {
        return ! ( ( *this < other ) || ( *this > other ) );
    }

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend bool operator== ( const Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( x == Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test on inequality
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if not equal, @c false otherwise
     */
    RATIONAL_CONSTEXPR inline bool operator!= ( const Rational &other ) const {
        return ! ( *this == other );
    }

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend bool operator!= ( const Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( x != Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test if less than
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if less than @c other, @c false otherwise
     */
    RATIONAL_CONSTEXPR bool operator< ( const Rational &other ) const;

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend bool operator< ( const Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( x < Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test if less or equal than
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if less or equal than @c other, @c false otherwise
     */
    RATIONAL_CONSTEXPR inline bool operator<= ( const Rational &other ) const {
        return ! ( other < *this );
    }

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend bool operator<= ( const Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( x <= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test if greater than
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if greater than @c other, @c false otherwise
     */
    RATIONAL_CONSTEXPR inline bool operator> ( const Rational &other ) const {
        return other < *this;
    }

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend bool operator> ( const Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( x > Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test if greater or equal than
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if greater or equal than @c other, @c false otherwise
     */
    RATIONAL_CONSTEXPR inline bool operator>= ( const Rational &other ) const {
        return ! ( *this < other );
    }

    template<template<typename, bool, template<class, typename, bool> class,
             template<typename> class> class U, template<class, typename, bool> class V>
    RATIONAL_CONSTEXPR inline friend bool operator>= ( const Rational& x,
            const Rational<integer_type, U, V>& y ) {
        return ( x >= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test if it is the neutral element to addition and subtraction
     *
     * Tests if the @c numerator is equal to the default constructed @c integer_type
     *
     * @note @f$ 0 + x = x + 0 = x @f$ and @f$ 0 - x = x - 0 = x @f$
     *
     * @warning it does @b not do the test for the neutral element to multiplication and division
     *
     * @return @c true if it is the neutral element to addition and subtraction, @c false otherwise
     */
    RATIONAL_CONSTEXPR inline bool operator!() const {
        return m_numer == zero_;
    }

    /**
     * @brief generates the string representation of %Rational
     *
     * @param[in] mixed if @c true, than a mixed (inproper) fraction is generated
     *
     * @return the string representation of %Rational
     */
    std::string str ( bool mixed = false ) const;

    /**
     * @brief output stream operator
     *
     * Sends a string representation of %Rational to a @c std::ostream
     *
     * @see str(bool)
     *
     * @param[out] o the stream to send to
     * @param[in] r the %Rational to send to the stream @c o
     *
     * @return the stream @c o
     */
    RATIONAL_CONSTEXPR friend std::ostream &operator<< ( std::ostream &o, const Rational &r ) {
        return ( o << r.str() );
    }

    /**
     * @brief reads in an @c expression from a @c std::istream and assign it to @c r
     *
     * @see Rational(const char *expr)
     *
     * @param[in] i the stream to read from
     * @param[out] r the %Rational to assign to
     *
     * @return the stream @c i
     */
    friend std::istream &operator>> ( std::istream &i, Rational &r ) {

        r = Rational ( std::string ( std::istreambuf_iterator<char> ( i ),
                                     std::istreambuf_iterator<char>() ).c_str() );
        return i;
    }

private:
    Rational &reduce();

    static std::size_t md ( integer_type &out, const typename std::vector<integer_type> &dv,
                            const integer_type &base );

    RATIONAL_CONSTEXPR inline static bool isLeftAssoc ( const char op ) {
        return op > 2;
    }

    RATIONAL_CONSTEXPR inline static unsigned char getPrec ( const char op ) {
        return !isLeftAssoc ( op ) ? 2 : ( ( op == '*' || op == '/' || op == '%' ) ? 1 : 0 );
    }

    typedef std::stack<Rational, std::vector<Rational> > evalStack;

    static bool eval ( const char op, evalStack &s, const char *expr );

    static void pushToken ( evalStack &rpn, const char *expr, std::ptrdiff_t tok_start,
                            std::ptrdiff_t *tok_len ) {

        rpn.push ( TYPE_CONVERT<const char *> ( expr + tok_start,
                                                expr + tok_start + *tok_len ).
                   template convert<typename
                   ExpressionEvalTraits<integer_type>::NumberType>() );

        *tok_len = 0;
    }

    template<class Op>
    Rational &knuth_addSub ( const Rational &o );

public:
    static const integer_type zero_; ///< represents @em zero in the given Rational::integer_type
    static const integer_type one_; ///< represents @em one in the given Rational::integer_type

private:
    integer_type m_numer;
    integer_type m_denom;
};

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
const typename Rational<T, GCD, CHKOP>::integer_type Rational<T, GCD, CHKOP>::zero_
    = integer_type();

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
const typename Rational<T, GCD, CHKOP>::integer_type Rational<T, GCD, CHKOP>::one_
    = integer_type ( 1 );

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP>::Rational ( const integer_type &n, const integer_type &d ) : m_numer ( n ),
    m_denom ( d ) {

#ifdef __EXCEPTIONS
    if ( m_denom == zero_ ) throw std::domain_error ( "denominator can't be null" );
#endif

    reduce();
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
template<typename NumberType> Rational<T, GCD, CHKOP>::Rational ( const NumberType &nt ) :
    m_numer ( TYPE_CONVERT<NumberType> ( nt ).template convert<integer_type> () ),
    m_denom ( one_ ) {

    _approxFract<integer_type, GCD, CHKOP, NumberType,
                 ! ( std::numeric_limits<NumberType>::is_integer ||
                     std::numeric_limits<NumberType>::is_exact ) >() ( *this, nt );
}

#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic push
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP>::Rational ( const rf_info &info ) : m_numer (), m_denom () {

    using namespace std;

    *this = ( Rational ( info.pre, info.reptend, info.reptend == zero_ ? one_ :
                         static_cast<integer_type> ( _type_round_helper<integer_type>() (
                                     pow10 ( ceil ( log10 ( ( info.reptend < integer_type() ?
                                             integer_type ( -info.reptend ) : info.reptend ) +
                                             one_ ) ) + info.leading_zeros ) - one_ ) ) ) *=
                  Rational ( one_, static_cast<integer_type> ( _type_round_helper<integer_type>() (
                                 pow10 ( ceil ( log10 ( ( info.pre < integer_type() ?
                                         integer_type ( -info.pre ) : info.pre ) + one_ ) ) +
                                         info.pre_leading_zeros ) ) ) ) );
}
#pragma GCC diagnostic pop

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP>::Rational ( const char *expr ) : m_numer(), m_denom ( one_ ) {

    if ( expr && *expr ) {

        static const char td_op[12] = { '\t', '\n', ' ', '(', ')', 1, 2, '%', '-', '/', '*', '+' };
        static const char *td_op_end = td_op + 12, *op_start = td_op + 5;

        std::stack<char, std::vector<char> > syard;
        std::ptrdiff_t tok_start = 0, tok_len = 0;
        evalStack rpn;

        const char *ptr = expr;
        char top, prev = 0;

        while ( *ptr ) {

            if ( std::find ( td_op, td_op_end, *ptr ) == td_op_end ) {

                if ( ( *ptr >= '0' && *ptr <= '9' ) || *ptr == '.' ) {
                    if ( !tok_len++ ) tok_start = ptr - expr;
                } else {
#ifdef __EXCEPTIONS
                    throw std::runtime_error ( std::string
                                               ( "invalid character(s) in expression: " ).
                                               append ( expr ) );
#endif
                }

                if ( ! * ( ptr + 1 ) ) pushToken ( rpn, expr, tok_start, &tok_len );

                prev = *ptr++;
                continue;

            } else if ( tok_len ) pushToken ( rpn, expr, tok_start, &tok_len );

            if ( *ptr == ' ' || *ptr == '\t' || *ptr == '\n' ) {
                ++ptr;
                continue;
            }

            if ( *ptr == '(' ) {

                prev = *ptr;
                syard.push ( *ptr );

            } else if ( *ptr == ')' ) {

                prev = *ptr;

                while ( !syard.empty() && ( top = syard.top() ) != '(' ) {

                    if ( !eval ( top, rpn, expr ) ) {
#ifdef __EXCEPTIONS
                        throw std::runtime_error ( std::string ( "invalid expression: " )
                                                   .append ( expr ) );
#endif
                    }

                    syard.pop();
                }

                if ( !syard.empty() && top == '(' ) {

                    syard.pop();

                } else {
#ifdef __EXCEPTIONS
                    throw std::runtime_error ( "mismatched braces" );
#endif
                }

            } else if ( std::find ( op_start, td_op_end, *ptr ) != td_op_end ) {

                char cop = *ptr;

                const bool isUnary = ( ptr == expr ) || ( prev == '(' ||
                                     std::find ( op_start, td_op_end, prev ) != td_op_end );

                if ( *ptr == '-' && isUnary ) {
                    cop = 1;
                } else if ( *ptr == '+' && isUnary ) {
                    cop = 2;
                } else {

                    while ( !syard.empty() &&
                            std::find ( op_start, td_op_end, top = syard.top() ) != td_op_end &&
                            ( ( isLeftAssoc ( cop ) && getPrec ( cop ) <= getPrec ( top ) ) ||
                              ( !isLeftAssoc ( cop ) && getPrec ( cop ) < getPrec ( top ) ) ) ) {

                        if ( !eval ( top, rpn, expr ) ) {
#ifdef __EXCEPTIONS
                            throw std::runtime_error ( std::string ( "invalid expression: " )
                                                       .append ( expr ) );
#endif
                        }

                        syard.pop();
                    }
                }

                prev = *ptr;
                syard.push ( cop );
            }

            ++ptr;
        }

        while ( !syard.empty() &&
                std::find ( op_start, td_op_end, top = syard.top() ) != td_op_end ) {

            if ( !eval ( top, rpn, expr ) ) {
#ifdef __EXCEPTIONS
                throw std::runtime_error ( std::string ( "invalid expression: " ).append ( expr ) );
#endif
            }

            syard.pop();
        }

        if ( ! ( !syard.empty() || rpn.empty() || rpn.size() > 1 ) ) {
            *this = rpn.top();
        } else {
#ifdef __EXCEPTIONS
            throw std::runtime_error ( std::string ( "invalid expression: " ).append ( expr ) );
#endif
        }
    }
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP>::~Rational() {}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP> &Rational<T, GCD, CHKOP>::reduce() {

    typename tmp::_ifThenElse<tmp::_isClassT<integer_type>::Yes, const integer_type &,
             const integer_type>::ResultT x ( m_numer != zero_ ?
                     GCD<T, std::numeric_limits<integer_type>::is_signed,
                     CHKOP, TYPE_CONVERT>() ( m_numer, m_denom ) : m_denom );

    if ( x != one_ ) {
        m_numer = op_divides() ( m_numer, x );
        m_denom = op_divides() ( m_denom, x );
    }

    return _swapSign<integer_type, GCD, CHKOP,
           std::numeric_limits<integer_type>::is_signed>() ( *this );
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
template<class Op> Rational<T, GCD, CHKOP> &
Rational<T, GCD, CHKOP>::knuth_addSub ( const Rational<T, GCD, CHKOP> &o ) {

    typename tmp::_ifThenElse<tmp::_isClassT<integer_type>::Yes, const integer_type &,
             const integer_type>::ResultT d1 ( GCD<integer_type,
                     std::numeric_limits<integer_type>::is_signed,
                     CHKOP, TYPE_CONVERT>() ( m_denom, o.m_denom ) );

    if ( d1 == one_ ) {

        m_numer = Op() ( op_multiplies() ( m_numer, o.m_denom ),
                         op_multiplies() ( m_denom, o.m_numer ) );
        m_denom = op_multiplies() ( m_denom, o.m_denom );

    } else {

        typename tmp::_ifThenElse<tmp::_isClassT<integer_type>::Yes, const integer_type &,
                 const integer_type>::ResultT t ( Op() (
                             op_multiplies() ( m_numer, ( op_divides() ( o.m_denom, d1 ) ) ),
                             op_multiplies() ( o.m_numer,
                                               ( op_divides() ( m_denom, d1 ) ) ) ) );

        typename tmp::_ifThenElse<tmp::_isClassT<integer_type>::Yes, const integer_type &,
                 const integer_type>::ResultT d2 ( GCD<integer_type,
                         std::numeric_limits<integer_type>::is_signed,
                         CHKOP, TYPE_CONVERT>() ( t, d1 ) );

        m_numer = op_divides() ( t, d2 );
        m_denom = op_multiplies() ( op_divides() ( m_denom, d1 ), op_divides() ( o.m_denom, d2 ) );
    }

    return *this;
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RATIONAL_CONSTEXPR Rational<T, GCD, CHKOP>
Rational<T, GCD, CHKOP>::operator+ ( const Rational& other ) const {
    return ( Rational ( *this ) += other );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> RATIONAL_CONSTEXPR inline
NumberType &operator+= ( NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( n = Rational<T, GCD, CHKOP> ( n ) += o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType> RATIONAL_CONSTEXPR inline
Rational<T, GCD, CHKOP> operator+ ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return ( o + Rational<T, GCD, CHKOP> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> RATIONAL_CONSTEXPR inline
Rational<T, GCD, CHKOP> operator+ ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( Rational<T, GCD, CHKOP> ( n ) + o );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> RATIONAL_CONSTEXPR inline
NumberType &operator-= ( NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( n = Rational<T, GCD, CHKOP> ( n ) -= o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType> RATIONAL_CONSTEXPR inline
Rational<T, GCD, CHKOP> operator- ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return o - Rational<T, GCD, CHKOP> ( n );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> RATIONAL_CONSTEXPR inline
Rational<T, GCD, CHKOP> operator- ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( Rational<T, GCD, CHKOP> ( n ) - o );
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RATIONAL_CONSTEXPR Rational<T, GCD, CHKOP>
Rational<T, GCD, CHKOP>::operator* ( const Rational& other ) const {
    return ( Rational ( *this ) *= other );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> RATIONAL_CONSTEXPR inline
NumberType &operator*= ( NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( n = Rational<T, GCD, CHKOP> ( n ) *= o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType> RATIONAL_CONSTEXPR inline
Rational<T, GCD, CHKOP> operator* ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return ( o * Rational<T, GCD, CHKOP> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> RATIONAL_CONSTEXPR inline
Rational<T, GCD, CHKOP> operator* ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( Rational<T, GCD, CHKOP> ( n ) * o );
}

template<typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP> &Rational<T, GCD, CHKOP>::invert() {

    using namespace std;
    swap ( m_numer, m_denom );

#ifdef __EXCEPTIONS
    if ( m_denom == zero_ ) throw std::domain_error ( "division by zero" );
#endif

    return _swapSign<integer_type, GCD, CHKOP,
           std::numeric_limits<integer_type>::is_signed>() ( *this );
}

template<typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> RATIONAL_CONSTEXPR
Rational<T, GCD, CHKOP> Rational<T, GCD, CHKOP>::inverse() const {
    return Rational ( *this ).invert();
}

template<typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> RATIONAL_CONSTEXPR
typename Rational<T, GCD, CHKOP>::mod_type Rational<T, GCD, CHKOP>::mod() const {
    return _mod<integer_type, GCD, CHKOP,
           std::numeric_limits<integer_type>::is_signed>() ( *this );
}

template<typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP>
std::size_t Rational<T, GCD, CHKOP>::md ( integer_type &out,
        const typename std::vector<integer_type> &dv, const integer_type &base ) {

    std::size_t zeros = 0u;

    if ( !dv.empty() ) {

        typename std::vector<integer_type>::const_iterator j ( dv.begin() );

        if ( ( out = *j++ ) == zero_ ) ++zeros;

        bool pre_zeros = false;

        for ( ; j != dv.end(); ++j ) {

            out = op_plus() ( op_multiplies() ( out, base ), *j );

            if ( !pre_zeros && *j == zero_ ) {
                ++zeros;
            } else {
                pre_zeros = true;
            }
        }
    }

    return zeros;
}

template<typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> struct _remquo {
    inline T operator() ( const T &x, const T &y, T &quo ) const {
        return typename Rational<T, GCD, CHKOP>::op_minus() ( x,
                typename Rational<T, GCD, CHKOP>::op_multiplies() ( y,
                        ( quo = typename Rational<T, GCD, CHKOP>::op_divides() ( x, y ) ) ) );
    }
};

template<template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> struct _remquo<int, GCD, CHKOP> {

    inline int operator() ( int x, int y, int &quo ) const {

        using namespace std;

        const div_t &d ( div ( x, y ) );

        quo = d.quot;
        return d.rem;
    }
};

template<template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> struct _remquo<long, GCD, CHKOP> {

    inline long operator() ( long x, long y, long &quo ) const {

        using namespace std;

        const ldiv_t &d ( ldiv ( x, y ) );

        quo = d.quot;
        return d.rem;
    }
};

#if (_XOPEN_SOURCE >= 600 || _ISOC99_SOURCE || _POSIX_C_SOURCE >= 200112L)
template<template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP>
struct _remquo<long long, GCD, CHKOP> {

    inline long long operator() ( long long x, long long y, long long &quo ) const {

        using namespace std;

        const lldiv_t &d ( lldiv ( x, y ) );

        quo = d.quot;
        return d.rem;
    }
};
#endif

#pragma GCC diagnostic ignored "-Wtype-limits"
#pragma GCC diagnostic ignored "-Wconversion"
#pragma GCC diagnostic push
template<typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP>
typename Rational<T, GCD, CHKOP>::integer_type
Rational<T, GCD, CHKOP>::decompose ( rf_info &rf_info, const integer_type &base ) const {

    using namespace std;

    std::vector<integer_type> rd, dg;
    integer_type d ( m_numer );

    do {

        integer_type q, aux;

        rd.push_back ( ( aux = _remquo<T, GCD, CHKOP> () ( d, m_denom, q ) ) < zero_ ?
                       integer_type ( -aux ) : aux );
        dg.push_back ( q < zero_ ? integer_type ( -q ) : q );

        d = op_multiplies() ( base, rd.back() );

    } while ( rd.back() != zero_ && !std::count ( rd.rbegin() + 1, rd.rend(), rd.back() ) );

    rf_info.reptend = rf_info.pre = zero_;
    rf_info.leading_zeros = rf_info.pre_leading_zeros = 0u;

    rf_info.reptend_digits.clear();
    rf_info.pre_digits.clear();

    const bool hasPer = rd.back() != zero_;
    const bool hasPre = !hasPer || rd.back() != rd.front();

    const typename std::vector<integer_type>::iterator &pivot ( dg.begin() + ( hasPre ?
            std::distance ( rd.begin(), std::find ( rd.begin(), rd.end(), rd.back() ) ) : 0 ) + 1 );

    if ( hasPre ) {
        rf_info.pre_digits.reserve ( static_cast<typename std::vector<integer_type>::size_type>
                                     ( std::distance ( dg.begin() + 1, pivot ) ) );
        rf_info.pre_digits.insert ( rf_info.pre_digits.end(), dg.begin() + 1, pivot );
        rf_info.pre_leading_zeros = md ( rf_info.pre, rf_info.pre_digits, base );
    }

    if ( hasPer ) {
        rf_info.reptend_digits.reserve ( static_cast<typename std::vector<integer_type>::size_type>
                                         ( std::distance ( pivot, dg.end() ) ) );
        rf_info.reptend_digits.insert ( rf_info.reptend_digits.end(), pivot, dg.end() );
        rf_info.leading_zeros = md ( rf_info.reptend, rf_info.reptend_digits, base );
    }

    const bool isNegative = m_numer < zero_;

    if ( isNegative ) {

        if ( !rf_info.pre_digits.empty() ) {
            rf_info.pre_digits.front() = integer_type ( -rf_info.pre_digits.front() );
        }

        if ( !rf_info.reptend_digits.empty() ) {
            rf_info.reptend_digits.front() = integer_type ( -rf_info.reptend_digits.front() );
        }
    }

    return isNegative ? integer_type ( -dg.front() ) : dg.front();
}
#pragma GCC diagnostic pop

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> RATIONAL_CONSTEXPR inline
NumberType &operator/= ( NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( n = Rational<T, GCD, CHKOP> ( n ) /= o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType> RATIONAL_CONSTEXPR inline
Rational<T, GCD, CHKOP> operator/ ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return ( o / Rational<T, GCD, CHKOP> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T, template<typename, bool,
         template<class, typename, bool> class, template<typename> class> class GCD,
         template<class, typename, bool> class CHKOP> RATIONAL_CONSTEXPR inline
Rational<T, GCD, CHKOP> operator/ ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( Rational<T, GCD, CHKOP> ( n ) / o );
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP>& Rational<T, GCD, CHKOP>::operator*= ( const Rational& other ) {

    typename tmp::_ifThenElse<tmp::_isClassT<integer_type>::Yes, const integer_type &,
             const integer_type>::ResultT d2 ( GCD<integer_type,
                     std::numeric_limits<integer_type>::is_signed, CHKOP, TYPE_CONVERT>()
                     ( m_denom, other.m_numer ) );

    typename tmp::_ifThenElse<tmp::_isClassT<integer_type>::Yes, const integer_type &,
             const integer_type>::ResultT d1 ( d2 == one_ ? GCD<integer_type,
                     std::numeric_limits<integer_type>::is_signed, CHKOP, TYPE_CONVERT>()
                     ( m_numer, other.m_denom ) : one_ );

    if ( ! ( d2 == one_ && d1 == one_ ) ) {

        m_numer = op_multiplies() ( ( op_divides() ( m_numer, d1 ) ),
                                    ( op_divides() ( other.m_numer, d2 ) ) );
        m_denom = op_multiplies() ( ( op_divides() ( m_denom, d2 ) ),
                                    ( op_divides() ( other.m_denom, d1 ) ) );

    } else {
        m_numer = op_multiplies() ( m_numer, other.m_numer );
        m_denom = op_multiplies() ( m_denom, other.m_denom );
    }

    return *this;
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP>& Rational<T, GCD, CHKOP>::operator%= ( const Rational& o ) {

    if ( m_denom != o.m_denom ) {

        typename tmp::_ifThenElse<tmp::_isClassT<integer_type>::Yes, const integer_type &,
                 const integer_type>::ResultT l ( _lcm<integer_type, GCD, CHKOP,
                         std::numeric_limits<integer_type>::is_signed>()
                         ( m_denom, o.m_denom ) );

        typename tmp::_ifThenElse<tmp::_isClassT<integer_type>::Yes, const integer_type &,
                 const integer_type>::ResultT a ( op_multiplies()
                         ( op_divides() ( l, o.m_denom ), o.m_numer ) );

        m_numer = op_modulus() ( op_plus() ( op_modulus() ( op_multiplies() ( op_divides()
                                             ( l, m_denom ), m_numer ), a ), a ), a );
        m_denom = l;

    } else {
        m_numer = op_modulus() ( op_plus() ( op_modulus() ( m_numer, o.m_numer ), o.m_numer ),
                                 o.m_numer );
    }

    return reduce();
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
std::string Rational<T, GCD, CHKOP>::str ( bool mixed ) const {

    std::ostringstream os;

    if ( mixed && m_denom != one_ ) {

        const mod_type &p ( mod() );

        if ( p.first != zero_ ) os << p.first << ' ';

        os << p.second.abs().str ( false );

    } else if ( mixed && m_denom == one_ ) {

        const mod_type &p ( mod() );

        os << ( p.first + p.second.abs().numerator() );

    } else {

        os << m_numer;

        if ( m_denom != one_ ) os << '/' << m_denom;
    }

    return os.str();
}

#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic push
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
bool Rational<T, GCD, CHKOP>::eval ( const char op, evalStack &s, const char *expr ) {

    if ( !s.empty() ) {

        Rational operand[2] = { s.top(), Rational() };

        s.pop();

        if ( op > 2 && s.empty() ) return false;

        switch ( op ) {
        case 1:
            s.push ( -operand[0] );
            return true;
        case 2:
            s.push ( operand[0] );
            return true;
        case '%':
            operand[1] = s.top();
            s.pop();
            s.push ( operand[1] %= operand[0] );
            return true;
        case '+':
            operand[1] = s.top();
            s.pop();
            s.push ( operand[1] += operand[0] );
            return true;
        case '-':
            operand[1] = s.top();
            s.pop();
            s.push ( operand[1] -= operand[0] );
            return true;
        case '*':
            operand[1] = s.top();
            s.pop();
            s.push ( operand[1] *= operand[0] );
            return true;
        case '/':
            if ( operand[0] == Rational() ) {
                throw std::domain_error ( std::string ( "division by zero in expression: " ).
                                          append ( expr ) );
            }
            operand[1] = s.top();
            s.pop();
            s.push ( operand[1] /= operand[0] );
            return true;
        }
    }

    return false;
}
#pragma GCC diagnostic pop

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T,
         template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RATIONAL_CONSTEXPR inline
NumberType &operator%= ( NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( n = Rational<T, GCD, CHKOP> ( n ) %= o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType> RATIONAL_CONSTEXPR inline
Rational<T, GCD, CHKOP> operator% ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return ( o % Rational<T, GCD, CHKOP> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T,
         template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RATIONAL_CONSTEXPR inline
Rational<T, GCD, CHKOP> operator% ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( Rational<T, GCD, CHKOP> ( n ) % o );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T,
         template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RATIONAL_CONSTEXPR inline
bool operator== ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( Rational<T, GCD, CHKOP> ( n ) == o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType> RATIONAL_CONSTEXPR inline
bool operator== ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return ( o == Rational<T, GCD, CHKOP> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T,
         template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RATIONAL_CONSTEXPR inline
bool operator!= ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ! ( Rational<T, GCD, CHKOP> ( n ) == o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType> RATIONAL_CONSTEXPR inline
bool operator!= ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return ! ( o == Rational<T, GCD, CHKOP> ( n ) );
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RATIONAL_CONSTEXPR bool Rational<T, GCD, CHKOP>::operator< ( const Rational &other ) const {
    return ( m_numer * other.m_denom ) < ( other.m_numer * m_denom );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T,
         template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RATIONAL_CONSTEXPR inline bool operator< ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( Rational<T, GCD, CHKOP> ( n ) < o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType>
RATIONAL_CONSTEXPR inline bool operator< ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return ( o < Rational<T, GCD, CHKOP> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T,
         template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RATIONAL_CONSTEXPR inline
bool operator<= ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ! ( o < Rational<T, GCD, CHKOP> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType> RATIONAL_CONSTEXPR inline
bool operator<= ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return ! ( Rational<T, GCD, CHKOP> ( n ) < o );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T,
         template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RATIONAL_CONSTEXPR inline bool operator> ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return o < Rational<T, GCD, CHKOP> ( n );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T,
         template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType> RATIONAL_CONSTEXPR inline
bool operator> ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return Rational<T, GCD, CHKOP> ( n ) < o;
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename NumberType, typename T,
         template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
RATIONAL_CONSTEXPR inline
bool operator>= ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ! ( Rational<T, GCD, CHKOP> ( n ) < o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType> RATIONAL_CONSTEXPR inline
bool operator>= ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return ! ( o < Rational<T, GCD, CHKOP> ( n ) );
}

template<typename NumberType, template<typename> class EPSILON>
struct _approxUtils {

    inline static bool approximated ( const NumberType &af, const NumberType &nt ) {
        return abs ( af - nt ) < eps_;
    }

    const static NumberType eps_;

private:
    inline static NumberType abs ( const NumberType &nt ) {
        return nt < NumberType() ? NumberType ( -nt ) : nt;
    }
};

template<typename NumberType, template<typename> class EPSILON>
const NumberType _approxUtils<NumberType, EPSILON>::eps_ ( EPSILON<NumberType>::value() );

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType, template<typename> class EPSILON, template<typename> class CONV>
struct _approxFract<T, GCD, CHKOP, NumberType, true, EPSILON, CONV> {
private:
    typedef Rational<T, GCD, CHKOP> rat;

public:
    void operator() ( rat &r, const NumberType &nt ) const;

private:
    const static T zero_;
    const static T one_;
};

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType, template<typename> class EPSILON, template<typename> class CONV>
const T _approxFract<T, GCD, CHKOP, NumberType, true, EPSILON, CONV>::zero_ = T();

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType, template<typename> class EPSILON, template<typename> class CONV>
const T _approxFract<T, GCD, CHKOP, NumberType, true, EPSILON, CONV>::one_ = T ( 1 );

#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic push
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType, template<typename> class EPSILON, template<typename> class CONV>
void _approxFract<T, GCD, CHKOP, NumberType, true, EPSILON, CONV>::operator() ( rat &r,
        const NumberType &nt ) const {

#ifdef __EXCEPTIONS
    if ( ! ( ! ( std::numeric_limits<T>::max() == zero_ ||
                 std::numeric_limits<T>::min() == zero_ ) &&
             ( nt > CONV<T> ( std::numeric_limits<T>::max() ).template convert<NumberType>() ||
               nt < CONV<T> ( std::numeric_limits<T>::min() ).template convert<NumberType>() ) ) ) {
#endif

        T m[2][2] = { { zero_, one_ }, { one_, zero_ } };

        NumberType x ( nt );

        typename tmp::_ifThenElse<tmp::_isClassT<T>::Yes, const NumberType &,
                 const NumberType>::ResultT one ( CONV<T> ( one_ ).template convert<NumberType>() );

        while ( !_approxUtils<NumberType, EPSILON>::approximated ( r, nt ) ) {

            using namespace std;

            typename tmp::_ifThenElse<tmp::_isClassT<T>::Yes, const T &,
                     const T>::ResultT n ( CONV<NumberType> ( floor ( x ) ).template convert<T>() );

            r.m_numer = typename rat::op_plus ()
                        ( m[0][0], typename rat::op_multiplies () ( n, m[0][1] ) );

            m[0][0] = m[0][1];
            m[0][1] = r.m_numer;

            r.m_denom = typename rat::op_plus ()
                        ( m[1][0], typename rat::op_multiplies () ( n, m[1][1] ) );

            m[1][0] = m[1][1];
            m[1][1] = r.m_denom;

            typename tmp::_ifThenElse<tmp::_isClassT<NumberType>::Yes, const NumberType &,
                     const NumberType>::ResultT
                     d ( x - CONV<T> ( n ).template convert<NumberType>() );

            if ( _approxUtils<NumberType, EPSILON>::approximated ( d, NumberType() ) ) break;

            x = one / d;
        }
#ifdef __EXCEPTIONS
    } else {
        throw std::domain_error ( "rational approximation overflow" );
    }
#endif
}
#pragma GCC diagnostic pop

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename NumberType, template<typename> class EPSILON, template<typename> class CONV>
struct _approxFract<T, GCD, CHKOP, NumberType, false, EPSILON, CONV> {
    inline void operator() ( const Rational<T, GCD, CHKOP> &,
                             const NumberType & ) const RATIONAL_NOEXCEPT {}
};

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _abs<T, GCD, CHKOP, true> {

    RATIONAL_CONSTEXPR inline Rational<T, GCD, CHKOP>
    operator() ( const Rational<T, GCD, CHKOP> &r ) const {
        return r.numerator() < Rational<T, GCD, CHKOP>::zero_ ? -r : r;
    }
};

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _abs<T, GCD, CHKOP, false> {

    RATIONAL_CONSTEXPR inline Rational<T, GCD, CHKOP>
    operator() ( const Rational<T, GCD, CHKOP> &r ) const {
        return r;
    }
};

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _mod<T, GCD, CHKOP, true> {

    typedef std::pair<T, Rational<T, GCD, CHKOP> > pair_type;

    pair_type operator() ( const Rational<T, GCD, CHKOP> &r ) const;
};

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
typename _mod<T, GCD, CHKOP, true>::pair_type
_mod<T, GCD, CHKOP, true>::operator() ( const Rational<T, GCD, CHKOP> &r ) const {

    const Rational<T, GCD, CHKOP> &h ( Rational<T, GCD, CHKOP> (
                                           ( typename Rational<T, GCD, CHKOP>::op_modulus()
                                                   ( r.m_numer, r.m_denom ) ), r.m_denom ) );

    return std::make_pair ( typename Rational<T, GCD, CHKOP>::op_divides() ( r.m_numer,
                            r.m_denom ), h );
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _mod<T, GCD, CHKOP, false> {

    typedef std::pair<T, Rational<T, GCD, CHKOP> > pair_type;

    RATIONAL_CONSTEXPR inline pair_type operator() ( const Rational<T, GCD, CHKOP> &r ) const {
        return std::make_pair ( typename Rational<T, GCD, CHKOP>::op_divides() ( r.m_numer,
                                r.m_denom ), Rational<T, GCD, CHKOP> (
                                    ( typename Rational<T, GCD, CHKOP>::op_modulus() ( r.m_numer,
                                            r.m_denom ) ), r.m_denom ) );
    }
};

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _pow<T, GCD, CHKOP, false> {

    inline Rational<T, GCD, CHKOP> operator() ( const Rational<T, GCD, CHKOP> &r,
            const T &exp ) const {
#ifdef __EXCEPTIONS
        if ( exp > Rational<T, GCD, CHKOP>::zero_ ) {
#endif

            Rational<T, GCD, CHKOP> b ( r );
            Rational<T, GCD, CHKOP> result ( Rational<T, GCD, CHKOP>::one_,
                                             Rational<T, GCD, CHKOP>::one_ );
            T e ( exp );

            do {

                if ( ( e & 1 ) != Rational<T, GCD, CHKOP>::zero_ ) result *= b;

                e >>= 1;
                b *= b;

            } while ( e != Rational<T, GCD, CHKOP>::zero_ );

            return result;

#ifdef __EXCEPTIONS
        } else {
            throw std::domain_error ( "power is undefined for zero" );
        }
#endif
    }
};

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _pow<T, GCD, CHKOP, true> {

    inline Rational<T, GCD, CHKOP> operator() ( const Rational<T, GCD, CHKOP> &r,
            const T &exp ) const {

#ifdef __EXCEPTIONS
        if ( exp >= Rational<T, GCD, CHKOP>::zero_ ) {
#endif
            return _pow<T, GCD, CHKOP, false>() ( r, exp );

#ifdef __EXCEPTIONS
        } else {
            throw std::domain_error ( "power is undefined for negative numbers" );
        }
#endif
    }
};

template<typename T, bool IsSigned, template<class, typename, bool> class CHKOP,
         template<typename> class CONV> struct GCD_null {

    RATIONAL_CONSTEXPR inline const T &operator() ( const T&, const T& ) const RATIONAL_NOEXCEPT {
        return one_;
    }

private:
    const static T one_;
};

template<typename T, bool IsSigned, template<class, typename, bool> class CHKOP,
         template<typename> class CONV>
const T GCD_null<T, IsSigned, CHKOP, CONV>::one_ = T ( 1 );

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_euclid_fast<T, false, CHKOP, CONV> {

    T operator() ( const T &a, const T &b ) const;

private:
    static const T zero_;
};

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
const T GCD_euclid_fast<T, false, CHKOP, CONV>::zero_ = T();

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
T GCD_euclid_fast<T, false, CHKOP, CONV>::operator() ( const T &a, const T &b ) const {

    T x ( a ), y ( b );

    while ( y != zero_ ) {

        x %= y;
        y ^= x;
        x ^= y;
        y ^= x;
    }

    return x;
}

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_euclid_fast<T, true, CHKOP, CONV> {

    T operator() ( const T &a, const T &b ) const;

private:
    static const T zero_;
};

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
const T GCD_euclid_fast<T, true, CHKOP, CONV>::zero_ = T();

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
T GCD_euclid_fast<T, true, CHKOP, CONV>::operator() ( const T &a, const T &b ) const {
    typename tmp::_ifThenElse<tmp::_isClassT<T>::Yes, const T &,
             const T>::ResultT h ( GCD_euclid_fast<T, false, CHKOP, CONV>() ( a, b ) );
    return h < zero_ ? T ( -h ) : h;
}

template<typename T, template<class, typename = T, bool = false> class CHKOP,
         template<typename> class CONV> struct GCD_euclid<T, false, CHKOP, CONV> {

    inline T operator() ( const T &a, const T &b ) const {

        T x ( a ), y ( b );

        while ( y != zero_ ) {

            typename tmp::_ifThenElse<tmp::_isClassT<T>::Yes, const T &,
                     const T>::ResultT h ( CHKOP<std::modulus<T> >() ( x, y ) );

            x = y;
            y = h;
        }

        return x;
    }

private:
    static const T zero_;
};

template<typename T, template<class, typename = T, bool = false> class CHKOP,
         template<typename> class CONV>
const T GCD_euclid<T, false, CHKOP, CONV>::zero_ = T();

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_euclid<T, true, CHKOP, CONV> {

    T operator() ( const T &a, const T &b ) const;

private:
    static const T zero_;
};

template<typename T, template<class, typename = T, bool = false> class CHKOP,
         template<typename> class CONV>
const T GCD_euclid<T, true, CHKOP, CONV>::zero_ = T();

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
T GCD_euclid<T, true, CHKOP, CONV>::operator() ( const T &a, const T &b ) const {
    typename tmp::_ifThenElse<tmp::_isClassT<T>::Yes, const T &,
             const T>::ResultT h ( GCD_euclid<T, false, CHKOP, CONV>() ( a, b ) );
    return h < zero_ ? T ( -h ) : h;
}

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_stein<T, false, CHKOP, CONV> {

    T operator() ( const T &a, const T &b ) const;

private:
    static const T zero_;
};

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
const T GCD_stein<T, false, CHKOP, CONV>::zero_ = T();

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
T GCD_stein<T, false, CHKOP, CONV>::operator() ( const T &a, const T &b ) const {

    T x ( a ), y ( b ), f = zero_;

    while ( y != zero_ ) {

        if ( x < y ) {

            y ^= x;
            x ^= y;
            y ^= x;

        } else if ( T ( x & 1 ) == zero_ ) {

            x >>= 1;

            if ( T ( y & 1 ) == zero_ ) {
                y >>= 1;
                ++f;
            }

        } else if ( T ( y & 1 ) == zero_ ) {
            y >>= 1;
        } else {
            x -= y;
        }
    }

    return x << CONV<T> ( f ).template convert<unsigned long int>();
}

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
struct GCD_stein<T, true, CHKOP, CONV> {

    T operator() ( const T &a, const T &b ) const;

private:
    static const T zero_;
};

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
const T GCD_stein<T, true, CHKOP, CONV>::zero_ = T();

template<typename T, template<class, typename, bool> class CHKOP, template<typename> class CONV>
T GCD_stein<T, true, CHKOP, CONV>::operator() ( const T &a, const T &b ) const {
    return GCD_stein<T, false, CHKOP, CONV>()
           ( a < zero_ ? T ( -a ) : a, b < zero_ ? T ( -b ) : b );
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _lcm<T, GCD, CHKOP, true> {
    T operator() ( const T &a, const T &b ) const;
};

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
T _lcm<T, GCD, CHKOP, true>::operator() ( const T &a, const T &b ) const {

    typename tmp::_ifThenElse<tmp::_isClassT<T>::Yes, const T &, const T>::ResultT
    x ( a < Rational<T, GCD, CHKOP>::zero_ ? T ( -a ) : a ),
      y ( b < Rational<T, GCD, CHKOP>::zero_ ? T ( -b ) : b );

    return typename Rational<T, GCD, CHKOP>::op_multiplies()
           ( ( typename Rational<T, GCD, CHKOP>::op_divides() ( static_cast<T> ( x ),
                   ( a != Rational<T, GCD, CHKOP>::zero_ ? GCD<T, false,
                     CHKOP, TYPE_CONVERT>() ( x, y ) : b ) ) ), static_cast<T> ( y ) );
}

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _lcm<T, GCD, CHKOP, false> {

    RATIONAL_CONSTEXPR inline T operator() ( const T &a, const T &b ) const {
        return typename Rational<T, GCD, CHKOP>::op_multiplies()
               ( typename Rational<T, GCD, CHKOP>::op_divides()
                 ( a, ( a ? GCD<T, false, CHKOP, TYPE_CONVERT>() ( a, b ) : b ) ), b );
    }
};

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _swapSign<T, GCD, CHKOP, true> {

    inline Rational<T, GCD, CHKOP> &operator() ( Rational<T, GCD, CHKOP> &r ) const {

        if ( r.m_denom < zero_ ) {
            r.m_numer = -r.m_numer;
            r.m_denom = -r.m_denom;
        }

        return r;
    }

private:
    const static T zero_;
};

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
const T _swapSign<T, GCD, CHKOP, true>::zero_ = T();

template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
struct _swapSign<T, GCD, CHKOP, false> {

    RATIONAL_CONSTEXPR inline Rational<T, GCD, CHKOP> &
    operator() ( Rational<T, GCD, CHKOP> &r ) const RATIONAL_NOEXCEPT {
        return r;
    }
};

#ifdef __EXCEPTIONS
template<typename T>
struct ENABLE_OVERFLOW_CHECK<std::plus<T>, T, true> {

    inline T operator() ( const T &x, const T& y ) const {

        if ( ! ( ( ( y > T() ) && ( x > ( std::numeric_limits<T>::max() - y ) ) ) ||
                 ( ( y < T() ) && ( x < ( std::numeric_limits<T>::min() - y ) ) ) ) ) {

            return std::plus<T>() ( x, y );
        }

        throw std::domain_error ( "addition overflow" );
    }
};
#endif

#ifdef __EXCEPTIONS
template<typename T>
struct ENABLE_OVERFLOW_CHECK<std::minus<T>, T, true> {

    inline T operator() ( const T &x, const T& y ) const {

        if ( ! ( ( y > T() && x < std::numeric_limits<T>::min() + y ) ||
                 ( y < T() && x > std::numeric_limits<T>::max() + y ) ) ) {

            return std::minus<T>() ( x, y );
        }

        throw std::domain_error ( "subtraction overflow" );
    }
};
#endif

#ifdef __EXCEPTIONS
template<typename T>
struct ENABLE_OVERFLOW_CHECK<std::negate<T>, T, true> {

    inline T operator() ( const T &x ) const {

        if ( x != std::numeric_limits<T>::min() ) {
            return std::negate<T>() ( x );
        }

        throw std::domain_error ( "negation overflow" );
    }
};
#endif

#ifdef __EXCEPTIONS
template<typename T>
struct ENABLE_OVERFLOW_CHECK<std::multiplies<T>, T, true> {
    T operator() ( const T &x, const T& y ) const;
};
#endif

#ifdef __EXCEPTIONS
template<typename T>
T ENABLE_OVERFLOW_CHECK<std::multiplies<T>, T, true>::operator() ( const T &x, const T& y ) const {

    bool overflow = false;

    if ( x > T() ) {
        if ( y > T() ) {
            if ( x > ( std::numeric_limits<T>::max() / y ) ) {
                overflow = true;
            }
        } else {
            if ( y < ( std::numeric_limits<T>::min() / x ) ) {
                overflow = true;
            }
        }
    } else {
        if ( y > T() ) {
            if ( x < ( std::numeric_limits<T>::min() / y ) ) {
                overflow = true;
            }
        } else {
            if ( ( x != T() ) && ( y < ( std::numeric_limits<T>::max() / x ) ) ) {
                overflow = true;
            }
        }
    }

    if ( overflow ) throw std::domain_error ( "multiplication overflow" );

    return std::multiplies<T>() ( x, y );
}
#endif

#ifdef __EXCEPTIONS
template<typename T>
struct ENABLE_OVERFLOW_CHECK<std::divides<T>, T, true> {
    T operator() ( const T &x, const T& y ) const;
};
#endif

#ifdef __EXCEPTIONS
template<typename T>
T ENABLE_OVERFLOW_CHECK<std::divides<T>, T, true>::operator() ( const T &x, const T& y ) const {

    if ( ! ( ( y == T() ) || ( ( x == std::numeric_limits<T>::min() ) && ( y == -1 ) ) ) ) {
        return std::divides<T>() ( x, y );
    }

    throw std::domain_error ( "division overflow" );
}
#endif

#ifdef __EXCEPTIONS
template<typename T>
struct ENABLE_OVERFLOW_CHECK<std::modulus<T>, T, true> {

    inline T operator() ( const T &x, const T& y ) const {

        if ( ! ( ( y == T() ) || ( ( x == std::numeric_limits<T>::min() ) && ( y == -1 ) ) ) ) {
            return std::modulus<T>() ( x, y );
        }

        throw std::domain_error ( "modulus overflow" );
    }
};
#endif

#ifdef __EXCEPTIONS
template<typename T>
struct ENABLE_OVERFLOW_CHECK<std::plus<T>, T, false> {

    inline T operator() ( const T &x, const T& y ) const {

        if ( ! ( std::numeric_limits<T>::max() - x < y ) ) {
            return std::plus<T>() ( x, y );
        }

        throw std::domain_error ( "unsigned addition wrap" );
    }
};
#endif

#ifdef __EXCEPTIONS
template<typename T>
struct ENABLE_OVERFLOW_CHECK<std::minus<T>, T, false> {

    inline T operator() ( const T &x, const T& y ) const {

        if ( ! ( x < y ) ) return std::minus<T>() ( x, y );

        throw std::domain_error ( "unsigned subtraction wrap" );
    }
};
#endif

template<typename T>
struct ENABLE_OVERFLOW_CHECK<std::negate<T>, T, false> {

#ifndef __EXCEPTIONS
    inline T operator() ( const T &x ) const {
        return std::negate<T>() ( x );
#else
    inline T operator() ( const T & ) const {
        throw std::domain_error ( "unsigned negation wrap" );
#endif
    }
};

#ifdef __EXCEPTIONS
template<typename T>
struct ENABLE_OVERFLOW_CHECK<std::multiplies<T>, T, false> {

    inline T operator() ( const T &x, const T& y ) const {

        if ( y == T ( 0 ) || ! ( x > std::numeric_limits<T>::max() / y ) ) {
            return std::multiplies<T>() ( x, y );
        }

        throw std::domain_error ( "unsigned multiplication wrap" );
    }
};
#endif

/**
 * @ingroup main
 * @brief Traits class to choose an appropriate Commons::Math::Rational
 */
template<typename T>
struct CFRationalTraits {
    typedef Rational<T> rational_type;
};

/**
 * @ingroup main
 * @brief Constructs a Rational of a given continued fraction sequence
 *
 * @see Commons::Math::CFRationalTraits
 *
 * @tparam IIter an input iterator
 *
 * @param[in] first iterator pointing to the begin of the sequence
 * @param[in] last iterator pointing to the end of the sequence
 *
 * @return the Commons::Math::Rational representing the sequence
 */
template<typename IIter>
typename CFRationalTraits<typename std::iterator_traits<IIter>::value_type>::rational_type
cf ( IIter first, IIter last ) {

    typedef typename std::iterator_traits<IIter>::value_type value_type;
    typedef typename CFRationalTraits<typename
    std::iterator_traits<IIter>::value_type>::rational_type rat;

    value_type m[2][2] = {
        { value_type(), value_type ( 1 ) },
        { value_type ( 1 ), value_type() }
    };

    value_type n = m[0][0], d = m[0][1];

    while ( first != last ) {

        n = typename rat::op_plus() (
                typename rat::op_multiplies() ( *first, m[0][1] ), m[0][0] );

        m[0][0] = m[0][1];
        m[0][1] = n;

        d = typename rat::op_plus() (
                typename rat::op_multiplies() ( *first++, m[1][1] ), m[1][0] );

        m[1][0] = m[1][1];
        m[1][1] = d;
    }

    return rat ( n, d );
}

/**
 * @ingroup main
 * @brief Extracts a continued fraction sequence of a Rational
 *
 * @see Commons::Math::Rational::mod()
 *
 * @tparam OIter an output iterator
 *
 * @param[in] r the Commons::Math::Rational to extract the continued fraction sequence from
 * @param[out] out iterator to ouput the sequence to
 *
 * @return iterator pointing to the end of @c out
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP,
         typename OIter> inline OIter seq ( const Rational<T, GCD, CHKOP> &r, OIter out ) {

    typedef Rational<T, GCD, CHKOP> rat;

    rat h ( r );

    typename rat::mod_type mt;

    do {

        mt = h.mod();
        * ( out++ ) = mt.first;

    } while ( mt.second.numerator() != T() && ( h = mt.second.invert(), true ) );

    return out;
}

}

}

namespace std {

/**
 * @ingroup main
 * @brief Overload of @c std::modf for @c %Rational types
 *
 * @see Rational::mod()
 *
 * @warning __iptr @b must point to a valid address and @b cannot be @c NULL
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 *
 * @param[in] __x the @c %Rational
 * @param[out] __iptr address of the @c integer_type to store the integral part
 *
 * @return the %Rational part of @c __x
 */
template<typename T, template<typename, bool, template<class, typename, bool> class,
         template<typename> class> class GCD, template<class, typename, bool> class CHKOP>
inline Commons::Math::Rational<T, GCD, CHKOP>
modf ( const Commons::Math::Rational<T, GCD, CHKOP> &__x,
       typename Commons::Math::Rational<T, GCD, CHKOP>::integer_type * __iptr ) {

    const typename Commons::Math::Rational<T, GCD, CHKOP>::mod_type &tmp ( __x.mod() );

    *__iptr = tmp.first;

    return tmp.second;
}

}

/**
 * @defgroup gcd Greatest common divisor algorithms
 *
 * The @em greatest @em common @em divisor @em algorithms are used to
 * reduce a Commons::Math::Rational so that
 * @f$ numerator \perp denominator @f$, i.e. @f$ gcd ( numerator, denominator ) = 1 @f$
 *
 * @b Example: \n A custom GCD algorithm could be implemented as following: @code{.cpp}
 * template<typename MyType, bool IsSigned, template<class, typename = MyType,
 *          bool = IsSigned> class CHKOP, template<typename = MyType> class CONV>
 *
 * struct GCD_myType {
 *
 *   inline MyType operator() ( const MyType &a, const MyType &b ) const {
 *       return myType_GCD_impl(a, b);
 *   }
 *
 * };@endcode and using it at second template parameter:
 * @code{.cpp} typedef Commons::Math::Rational<MyType, GCD_myType> myType_rational;@endcode
 */

#endif /* COMMONS_MATH_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
