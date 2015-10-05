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

#ifndef COMMONS_MATH_RATIONAL_H
#define COMMONS_MATH_RATIONAL_H

#include <functional>
#include <sstream>
#include <limits>
#include <cmath>

#ifdef __EXCEPTIONS
#include <stdexcept>
#endif

namespace Commons {

namespace Math {

template<typename, template<typename, bool, template<class, typename, bool> class> class,
         template<class, typename, bool> class, bool> struct _mod;

template<typename, template<typename, bool, template<class, typename, bool> class> class,
         template<class, typename, bool> class, bool> struct _abs;

template<typename, template<typename, bool, template<class, typename, bool> class> class,
         template<class, typename, bool> class, bool> struct _lcm;

template<template<typename, bool, template<class, typename, bool> class> class,
         template<class, typename, bool> class, bool> struct _changeSign;

template<typename T> struct TYPE_CONVERT {
    inline explicit TYPE_CONVERT ( const T& v ) : val ( v ) {}

    template<typename U> inline U convert() const {
        return static_cast<U> ( val );
    }

private:
    const T& val;
};

template<typename T> struct EPSILON {
    inline static T value() {
        return std::numeric_limits<T>::epsilon();
    }
};

template<typename, template<typename, bool, template<class, typename, bool> class> class,
         template<class, typename, bool> class, typename, bool,
         template<typename> class = EPSILON,
         template<typename> class = TYPE_CONVERT> struct _approxFract;

/**
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

    inline T operator() ( const T &x, const T& y ) const {
        return Op() ( x, y );
    }
};

template<typename T, bool IsSigned>
struct NO_OPERATOR_CHECK<std::negate<T>, T, IsSigned> {

    inline T operator() ( const T &x ) const {
        return std::negate<T>() ( x );
    }
};

/**
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

    inline T operator() ( const T &x, const T& y ) const {
        return NO_OPERATOR_CHECK<Op, T, IsSigned>() ( x, y );
    }
};

/**
 * @brief Stein GCD algorithm implementation
 *
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, bool IsSigned, template<class, typename = T, bool = IsSigned> class CHKOP>
struct GCD_stein;

/**
 * @brief Euclid GCD algorithm (safe) implementation
 *
 * This implementation supports overlow/wrap checking
 *
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, bool IsSigned, template<class, typename = T, bool = IsSigned> class CHKOP>
struct GCD_euclid;

/**
 * @brief Euclid GCD algorithm (fast) implementation
 *
 * @see GCD_euclid if your number class doesn't support all needed operators
 *
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, bool IsSigned, template<class, typename = T, bool = IsSigned> class CHKOP>
struct GCD_euclid_fast;

/**
 * @brief %Rational (fraction) template class
 *
 * @note All %Rational objects are reduced
 *
 * @tparam T storage type
 * @tparam GCD GCD algorithm
 * @tparam CHKOP checked operator @see ENABLE_OVERFLOW_CHECK
 */
template<typename T, template<typename, bool, template<class, typename, bool> class>
class GCD = GCD_euclid_fast,
      template<class, typename = T, bool = std::numeric_limits<T>::is_signed>
      class CHKOP = NO_OPERATOR_CHECK>
class Rational {
    friend struct _changeSign<GCD, CHKOP, std::numeric_limits<T>::is_signed>;
    friend struct _mod<T, GCD, CHKOP, std::numeric_limits<T>::is_signed>;
    template<typename, template<typename, bool, template<class, typename, bool> class> class,
             template<class, typename, bool> class, typename, bool, template<typename> class,
             template<typename> class> friend struct _approxFract;
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
     */
    Rational() : m_numer (), m_denom ( 1 ) {}

    /**
     * @brief copy constructor
     *
     * @param[in] other the %Rational to copy
     */
    Rational ( const Rational &other ) : m_numer ( other.m_numer ), m_denom ( other.m_denom ) {}

    /**
    * @brief creates a %Rational
    *
    * Creates a copy of @c numer and divides it by @c denom
    *
    * @note to achieve sth like 1/(1/2) you must explicitely cast the numerator, i.e.
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
    template<template<typename, bool, template<class, typename, bool> class> class U1,
             template<class, typename, bool> class V1,
             template<typename, bool, template<class, typename, bool> class> class U2,
             template<class, typename, bool> class V2>
    Rational ( const Rational<integer_type, U1, V1> &numer,
               const Rational<integer_type, U2, V2> &denom ) :
        m_numer ( numer.numerator() ), m_denom ( numer.denominator() ) {
        *this *= denom.inverse();
    }

    /**
     * @brief creates a %Rational
     *
     * @param[in] numer the numerator
     * @param[in] denom the denominator
     */
    Rational ( const integer_type &numer, const integer_type &denom );

    /**
     * @brief creates a inproper (mixed) %Rational
     *
     * @param[in] whole whole number part
     * @param[in] numer the numerator
     * @param[in] denom the denominator
     */
    Rational ( const integer_type &whole, const integer_type &numer,
               const integer_type &denom ) : m_numer (), m_denom ( 1 ) {
        *this += Rational ( numer, denom ) += Rational ( whole );
    }

    /**
     * @brief creates an approximated %Rational
     *
     * @tparam NumberType type of the number to approximate
     *
     * @param[in] number the number to create an approximated %Rational of
     */
    template<typename NumberType>
    Rational ( const NumberType &number );

    ~Rational();

    /**
     * @brief assign another %Rational
     *
     * @param[in] another the %Rational to assign
     */
    Rational &operator= ( const Rational& another ) {

        if ( this != &another ) {
            m_numer = another.m_numer;
            m_denom = another.m_denom;
        }

        return *this;
    }

    /**
     * @brief assigns from a @c NumberType
     *
     * The Number is approximated to a %Rational, then it gets assigned
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
    inline operator void() const {}
#endif

    /**
     * @brief convert to @c NumberType
     *
     * @tparam NumberType type of the number to approximate
     *
     * @return the number value of the %Rational
     */
    template<typename NumberType>
    inline operator NumberType() const {
        return static_cast<NumberType> ( m_numer ) / static_cast<NumberType> ( m_denom );
    }

    /**
     * @brief gets the numerator
     *
     * @return the numerator
     */
    inline integer_type numerator() const throw() {
        return m_numer;
    }

    /**
     * @brief gets the denominator
     *
     * @return the denominator
     */
    inline integer_type denominator() const throw() {
        return m_denom;
    }

    /**
     * @brief extract the integral and fractional part
     *
     * @return a @c mod_type containing the integral and fractional part
     */
    inline mod_type mod() const {
        return _mod<integer_type, GCD, CHKOP,
               std::numeric_limits<integer_type>::is_signed>() ( *this );
    }

    /**
     * @brief gets the absolute %Rational
     *
     * @return a copy of the absolute %Rational
     */
    inline Rational abs() const {
        return _abs<integer_type, GCD, CHKOP,
               std::numeric_limits<integer_type>::is_signed>() ( *this );
    }

    /**
     * @brief inverts the %Rational
     *
     * @return the inverted %Rational
     */
    Rational &invert() {

        using namespace std;
        swap ( m_numer, m_denom );

#ifdef __EXCEPTIONS
        if ( m_denom == integer_type() ) throw std::domain_error ( "division by zero" );
#endif

        return *this;
    }

    /**
     * @brief gets a copy of the inverted %Rational
     *
     * @return a copy of the inverted %Rational
     */
    inline Rational inverse() const {
        return Rational ( *this ).invert();
    }

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

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend Rational &operator+= ( Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( x += Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief add a %Rational
     *
     * @param[in] other the %Rational to add
     *
     * @return a new %Rational
     */
    Rational operator+ ( const Rational& other ) const;

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend Rational operator+ ( const Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( Rational ( x ) += Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief get a copy of the %Rational
     *
     * @return a copy of %Rational
     */
    inline Rational operator+() const {
        return Rational ( *this );
    }

    /**
     * @brief pre-increment a %Rational
     *
     * the result will be @code (numerator + denominator) / denominator @endcode
     *
     * @return the incremented %Rational
     */
    inline Rational& operator++() {
        m_numer = op_plus() ( m_numer, m_denom );
        return gcm ( *this );
    }

    /**
     * @brief post-increment a %Rational
     *
     * the result will be @code (numerator + denominator) / denominator @endcode
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

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend Rational &operator-= ( Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( x -= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief subtract a %Rational
     *
     * @param[in] other the %Rational to subtract
     *
     * @return a new %Rational
     */
    inline Rational operator- ( const Rational& other ) const {
        return ( Rational ( *this ) -= other );
    }

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend Rational operator- ( const Rational& x, const Rational<integer_type, U, V>& y ) {
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
     * the result will be @code (numerator - denominator) / denominator @endcode
     *
     * @return the decremented %Rational
     */
    inline Rational& operator--() {
        m_numer = op_minus() ( m_numer, m_denom );
        return gcm ( *this );
    }

    /**
     * @brief post-decrement a %Rational
     *
     * the result will be @code (numerator - denominator) / denominator @endcode
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

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend Rational &operator*= ( Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( x *= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief multiply a %Rational
     *
     * @param[in] other the %Rational to multiply
     *
     * @return a new %Rational
     */
    Rational operator* ( const Rational& other ) const;

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend Rational operator* ( const Rational& x, const Rational<integer_type, U, V>& y ) {
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

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend Rational &operator/= ( Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( x /= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief divide a %Rational
     *
     * @param[in] other the %Rational to divide
     *
     * @return a new %Rational
     */
    inline Rational operator/ ( const Rational& other ) const {
        return ( Rational ( *this ) /= other );
    }

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend Rational operator/ ( const Rational& x, const Rational<integer_type, U, V>& y ) {
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

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend Rational &operator%= ( Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( x %= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief modulo a %Rational
     *
     * @param[in] other the %Rational to modulo
     *
     * @return a new %Rational
     */
    inline Rational operator% ( const Rational& other ) const {
        return ( Rational ( *this ) %= other );
    }

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend Rational operator% ( const Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( Rational ( x ) %= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test on equality
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if equal, @c false otherwise
     */
    inline bool operator== ( const Rational &other ) const {
        return ! ( ( *this < other ) || ( *this > other ) );
    }

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend bool operator== ( const Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( x == Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test on inequality
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if not equal, @c false otherwise
     */
    inline bool operator!= ( const Rational &other ) const {
        return ! ( *this == other );
    }

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend bool operator!= ( const Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( x != Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test if less than
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if less than @c other, @c false otherwise
     */
    inline bool operator< ( const Rational &other ) const {
        return /* ( m_denom * o.m_denom ) > 0 ? */ ( m_numer * other.m_denom ) <
                ( other.m_numer * m_denom )
                /* : ( o.m_numer * m_denom ) < ( m_numer * o.m_denom ) */;
        // denom can NEVER be zero!
    }

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend bool operator< ( const Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( x < Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test if less or equal than
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if less or equal than @c other, @c false otherwise
     */
    inline bool operator<= ( const Rational &other ) const {
        return ! ( other < *this );
    }

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend bool operator<= ( const Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( x <= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test if greater than
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if greater than @c other, @c false otherwise
     */
    inline bool operator> ( const Rational &other ) const {
        return other < *this;
    }

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend bool operator> ( const Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( x > Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test if greater or equal than
     *
     * @param[in] other the %Rational to test to
     *
     * @return @c true if greater or equal than @c other, @c false otherwise
     */
    inline bool operator>= ( const Rational &other ) const {
        return ! ( *this < other );
    }

    template<template<typename, bool, template<class, typename, bool> class> class U,
             template<class, typename, bool> class V>
    inline friend bool operator>= ( const Rational& x, const Rational<integer_type, U, V>& y ) {
        return ( x >= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief test if it is the neutral element to addition and subtraction
     *
     * Tests if the @c numerator is equal to the default constructed @c integer_type
     *
     * @return @c true if it is the neutral element to addition and subtraction, @c false otherwise
     */
    inline bool operator!() const {
        return m_numer == integer_type();
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
    friend std::ostream &operator<< ( std::ostream &o, const Rational &r ) {
        return ( o << r.str() );
    }

    /**
     * @brief reads in a @c double from a @c std::istream and assign it to @c r
     *
     * @param[in] i the stream to read from
     * @param[out] r the %Rational to assign to
     *
     * @return the stream @c i
     */
    friend std::istream &operator>> ( std::istream &i, Rational &r ) {

        double d;

        i >> d;
        r = d;

        return i;
    }

private:
    Rational &gcm ( const Rational &o );

    template<class Op>
    Rational &knuth_addSub ( const Rational &o );

private:
    integer_type m_numer;
    integer_type m_denom;
};

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP>::Rational ( const integer_type &n, const integer_type &d ) : m_numer ( n ),
    m_denom ( d ) {

#ifdef __EXCEPTIONS
    if ( m_denom == integer_type() ) throw std::domain_error ( "denominator can't be null" );
#endif

    gcm ( *this );
}

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
template<typename NumberType> Rational<T, GCD, CHKOP>::Rational ( const NumberType &nt ) :
    m_numer ( static_cast<integer_type> ( nt ) ), m_denom ( 1 ) {

    _approxFract<integer_type, GCD, CHKOP, NumberType,
                 ! ( std::numeric_limits<NumberType>::is_integer ||
                     std::numeric_limits<NumberType>::is_exact ) >() ( *this, nt );
}

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP>::~Rational() {}

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP> &Rational<T, GCD, CHKOP>::gcm ( const Rational &o ) {

    const integer_type &x ( o.m_numer != integer_type() ?
                            GCD<T, std::numeric_limits<integer_type>::is_signed,
                            CHKOP>() ( o.m_numer, o.m_denom ) : o.m_denom );

    if ( x != static_cast<integer_type> ( 1 ) ) {
        m_numer = op_divides() ( m_numer, x );
        m_denom = op_divides() ( m_denom, x );
    }

    return _changeSign<GCD, CHKOP, std::numeric_limits<integer_type>::is_signed>() ( *this );
}

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP> template<class Op> Rational<T, GCD, CHKOP> &
Rational<T, GCD, CHKOP>::knuth_addSub ( const Rational<T, GCD, CHKOP> &o ) {

    const integer_type &d1 ( GCD<integer_type, std::numeric_limits<integer_type>::is_signed,
                             CHKOP>() ( m_denom, o.m_denom ) );

    if ( d1 == static_cast<integer_type> ( 1 ) ) {

        m_numer = Op() ( op_multiplies() ( m_numer, o.m_denom ),
                         op_multiplies() ( m_denom, o.m_numer ) );
        m_denom = op_multiplies() ( m_denom, o.m_denom );

    } else {

        const integer_type &t ( Op() (
                                    op_multiplies() ( m_numer, ( op_divides() ( o.m_denom, d1 ) ) ),
                                    op_multiplies() ( o.m_numer,
                                            ( op_divides() ( m_denom, d1 ) ) ) ) );

        const integer_type &d2 ( GCD<integer_type, std::numeric_limits<integer_type>::is_signed,
                                 CHKOP>() ( t, d1 ) );

        m_numer = op_divides() ( t, d2 );
        m_denom = op_multiplies() ( op_divides() ( m_denom, d1 ), op_divides() ( o.m_denom, d2 ) );
    }

    return *this;
}

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP> Rational<T, GCD, CHKOP>::operator+ ( const Rational& other ) const {
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
         template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline NumberType &operator+= ( NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType>
inline Rational<T, GCD, CHKOP> operator+ ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
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
         template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline Rational<T, GCD, CHKOP> operator+ ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
         template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline NumberType &operator-= ( NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType>
inline Rational<T, GCD, CHKOP> operator- ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
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
         template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline Rational<T, GCD, CHKOP> operator- ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( Rational<T, GCD, CHKOP> ( n ) - o );
}

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP> Rational<T, GCD, CHKOP>::operator* ( const Rational& other ) const {
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
         template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline NumberType &operator*= ( NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType>
inline Rational<T, GCD, CHKOP> operator* ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
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
         template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline Rational<T, GCD, CHKOP> operator* ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( Rational<T, GCD, CHKOP> ( n ) * o );
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
         template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline NumberType &operator/= ( NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType>
inline Rational<T, GCD, CHKOP> operator/ ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
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
         template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline Rational<T, GCD, CHKOP> operator/ ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
    return ( Rational<T, GCD, CHKOP> ( n ) / o );
}

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP>& Rational<T, GCD, CHKOP>::operator*= ( const Rational& other ) {

    const integer_type &d1 ( GCD<integer_type,
                             std::numeric_limits<integer_type>::is_signed, CHKOP>() ( m_numer,
                                     other.m_denom ) );
    const integer_type &d2 ( GCD<integer_type,
                             std::numeric_limits<integer_type>::is_signed, CHKOP>() ( m_denom,
                                     other.m_numer ) );

    if ( ! ( d1 == static_cast<integer_type> ( 1 ) &&
             d2 == static_cast<integer_type> ( 1 ) ) ) {

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

template<typename T,
         template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
Rational<T, GCD, CHKOP>& Rational<T, GCD, CHKOP>::operator%= ( const Rational& o ) {

    if ( m_denom != o.m_denom ) {

        const integer_type &l ( _lcm<integer_type, GCD, CHKOP,
                                std::numeric_limits<integer_type>::is_signed>()
                                ( m_denom, o.m_denom ) );

        const integer_type &a ( op_multiplies() ( op_divides() ( l, o.m_denom ), o.m_numer ) );

        m_numer = op_modulus() ( op_plus() ( op_modulus() ( op_multiplies() ( op_divides()
                                             ( l, m_denom ), m_numer ), a ), a ), a );
        m_denom = l;

    } else {
        m_numer = op_modulus() ( op_plus() ( op_modulus() ( m_numer, o.m_numer ), o.m_numer ),
                                 o.m_numer );
    }

    return gcm ( *this );
}

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
std::string Rational<T, GCD, CHKOP>::str ( bool mixed ) const {

    std::ostringstream os;

    if ( mixed && m_denom != T ( 1 ) ) {

        const mod_type &p ( mod() );

        if ( p.first != integer_type() ) os << p.first << ' ';

        os << p.second.str ( false );

    } else if ( mixed && m_denom == T ( 1 ) ) {

        const mod_type &p ( mod() );

        os << ( p.first + p.second.numerator() );

    } else {

        os << m_numer;

        if ( m_denom != T ( 1 ) ) os << '/' << m_denom;
    }

    return os.str();
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
         template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline NumberType &operator%= ( NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType>
inline Rational<T, GCD, CHKOP> operator% ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
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
         template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline Rational<T, GCD, CHKOP> operator% ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
         template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline bool operator== ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType>
inline bool operator== ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
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
         template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline bool operator!= ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
template<typename T,
         template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType>
inline bool operator!= ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return ! ( o == Rational<T, GCD, CHKOP> ( n ) );
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
         template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline bool operator< ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType>
inline bool operator< ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
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
         template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline bool operator<= ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType>
inline bool operator<= ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
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
         template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline bool operator> ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
         template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType>
inline bool operator> ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
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
         template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
inline bool operator>= ( const NumberType &n, const Rational<T, GCD, CHKOP>& o ) {
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
template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType>
inline bool operator>= ( const Rational<T, GCD, CHKOP>& o, const NumberType &n ) {
    return ! ( o < Rational<T, GCD, CHKOP> ( n ) );
}

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType,
         template<typename> class EPSILON, template<typename> class CONV>
struct _approxFract<T, GCD, CHKOP, NumberType, true, EPSILON, CONV> {
private:
    typedef Rational<T, GCD, CHKOP> rat;

public:
    void operator() ( rat &r, const NumberType &nt ) const;

private:
    inline NumberType abs ( const NumberType &nt ) const {
        return nt < NumberType() ? NumberType ( -nt ) : nt;
    }
};

#pragma GCC diagnostic ignored "-Wfloat-equal"
#pragma GCC diagnostic push
template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType,
         template<typename> class EPSILON, template<typename> class CONV>
void _approxFract<T, GCD, CHKOP, NumberType, true, EPSILON, CONV>::operator() ( rat &r,
        const NumberType &nt ) const {

#ifdef __EXCEPTIONS
    if ( ! ( ! ( std::numeric_limits<T>::max() == T() || std::numeric_limits<T>::min() == T() ) &&
             ( nt > std::numeric_limits<T>::max() || nt < std::numeric_limits<T>::min() ) ) ) {
#endif

        const NumberType &eps ( EPSILON<NumberType>::value() );

        T p[2] = { T(), T ( 1 ) };
        T q[2] = { T ( 1 ), T() };

        NumberType x ( nt );

        while ( ! ( abs ( ( CONV<T> ( r.m_numer ).template convert<NumberType>() /
                            CONV<T> ( r.m_denom ).template convert<NumberType>() ) - nt ) < eps ) ) {

            const T &n ( static_cast<T> ( std::floor ( x ) ) );
            x = static_cast<NumberType> ( 1 ) / ( x -  CONV<T> ( n ).template convert<NumberType>() );

            r.m_numer = typename rat::op_plus ()
                        ( p[0], typename rat::op_multiplies () ( n, p[1] ) );
            p[0] = p[1];
            p[1] = r.m_numer;

            r.m_denom = typename rat::op_plus ()
                        ( q[0], typename rat::op_multiplies () ( n, q[1] ) );
            q[0] = q[1];
            q[1] = r.m_denom;
        }
#ifdef __EXCEPTIONS
    } else {
        throw std::domain_error ( "rational approximation overflow" );
    }
#endif
}
#pragma GCC diagnostic pop

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP, typename NumberType,
         template<typename> class EPSILON, template<typename> class CONV>
struct _approxFract<T, GCD, CHKOP, NumberType, false, EPSILON, CONV> {
    inline void operator() ( const Rational<T, GCD, CHKOP> &, const NumberType & ) const {}
};

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP> struct _abs<T, GCD, CHKOP, true> {

    inline Rational<T, GCD, CHKOP> operator() ( const Rational<T, GCD, CHKOP> &r ) const {
        return r.numerator() < T() ? -r : r;
    }
};

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP> struct _abs<T, GCD, CHKOP, false> {

    inline Rational<T, GCD, CHKOP> operator() ( const Rational<T, GCD, CHKOP> &r ) const {
        return r;
    }
};

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP> struct _mod<T, GCD, CHKOP, true> {

    typedef std::pair<T, Rational<T, GCD, CHKOP> > pair_type;

    inline pair_type operator() ( const Rational<T, GCD, CHKOP> &r ) const {

        const Rational<T, GCD, CHKOP> &h ( Rational<T, GCD, CHKOP> (
                                               ( typename Rational<T, GCD, CHKOP>::op_modulus()
                                                       ( r.m_numer, r.m_denom ) ), r.m_denom ) );

        return std::make_pair ( typename Rational<T, GCD, CHKOP>::op_divides() ( r.m_numer,
                                r.m_denom ), r.m_numer < T() ? -h : h );
    }
};

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP> struct _mod<T, GCD, CHKOP, false> {

    typedef std::pair<T, Rational<T, GCD, CHKOP> > pair_type;

    inline pair_type operator() ( const Rational<T, GCD, CHKOP> &r ) const {
        return std::make_pair ( typename Rational<T, GCD, CHKOP>::op_divides() ( r.m_numer,
                                r.m_denom ), Rational<T, GCD, CHKOP> (
                                    ( typename Rational<T, GCD, CHKOP>::op_modulus() ( r.m_numer,
                                            r.m_denom ) ), r.m_denom ) );
    }
};

template<typename T, template<class, typename, bool> class CHKOP>
struct GCD_euclid_fast<T, false, CHKOP> {

    inline T operator() ( const T &a, const T &b ) const {

        T x ( a ), y ( b );

        // while ( y ) { const integer_type &h ( x % y ); x = y; y = h; }

        while ( y != T() ) {

            x %= y;
            y ^= x;
            x ^= y;
            y ^= x;
        }

        return x;
    }
};

template<typename T, template<class, typename, bool> class CHKOP>
struct GCD_euclid_fast<T, true, CHKOP> {
    T operator() ( const T &a, const T &b ) const;
};

template<typename T, template<class, typename, bool> class CHKOP>
T GCD_euclid_fast<T, true, CHKOP>::operator() ( const T &a, const T &b ) const {
    const T &h ( GCD_euclid_fast<T, false, CHKOP>() ( a, b ) );
    return h < T() ? T ( -h ) : h;
}

template<typename T, template<class, typename = T, bool = false> class CHKOP>
struct GCD_euclid<T, false, CHKOP> {

    inline T operator() ( const T &a, const T &b ) const {

        T x ( a ), y ( b );

        while ( y != T() ) {

            const T &h ( CHKOP<std::modulus<T> >() ( x, y ) );
            x = y;
            y = h;
        }

        return x;
    }
};

template<typename T, template<class, typename, bool> class CHKOP>
struct GCD_euclid<T, true, CHKOP> {
    T operator() ( const T &a, const T &b ) const;
};

template<typename T, template<class, typename, bool> class CHKOP>
T GCD_euclid<T, true, CHKOP>::operator() ( const T &a, const T &b ) const {
    const T &h ( GCD_euclid<T, false, CHKOP>() ( a, b ) );
    return h < T() ? T ( -h ) : h;
}

template<typename T, template<class, typename, bool> class CHKOP>
struct GCD_stein<T, false, CHKOP> {

    inline T operator() ( const T &a, const T &b ) const {

        T x ( a ), y ( b ), f = T();

        while ( y ) {

            if ( x < y ) {

                y ^= x;
                x ^= y;
                y ^= x;

            } else if ( ! ( x & 1 ) ) {

                x >>= 1;

                if ( ! ( y & 1 ) ) {
                    y >>= 1;
                    ++f;
                }

            } else if ( ! ( y & 1 ) ) {
                y >>= 1;
            } else {
                x -= y;
            }
        }

        return x << f;
    }
};

template<typename T, template<class, typename, bool> class CHKOP>
struct GCD_stein<T, true, CHKOP> {

    inline T operator() ( const T &a, const T &b ) const {
        return GCD_stein<T, false, CHKOP>() ( a < T() ? T ( -a ) : a, b < T() ? T ( -b ) : b );
    }
};

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP> struct _lcm<T, GCD, CHKOP, true> {

    inline T operator() ( const T &a, const T &b ) const {

        const T &x ( a < T() ? -a : a ), &y ( b < T() ? -b : b );

        return typename Rational<T, GCD, CHKOP>::op_multiplies()
               ( ( typename Rational<T, GCD, CHKOP>::op_divides() ( static_cast<T> ( x ),
                       ( a ? GCD<T, false, CHKOP>() ( x, y ) : b ) ) ), static_cast<T> ( y ) );
    }
};

template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP> struct _lcm<T, GCD, CHKOP, false> {

    inline T operator() ( const T &a, const T &b ) const {
        return typename Rational<T, GCD, CHKOP>::op_multiplies()
               ( typename Rational<T, GCD, CHKOP>::op_divides()
                 ( a, ( a ? GCD<T, false, CHKOP>() ( a, b ) : b ) ), b );
    }
};

template<template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
struct _changeSign<GCD, CHKOP, true> {

    template<typename T>
    inline Rational<T, GCD, CHKOP> &operator() ( Rational<T, GCD, CHKOP> &r ) const {

        if ( r.m_denom < T() ) {
            r.m_numer = -r.m_numer;
            r.m_denom = -r.m_denom;
        }

        return r;
    };
};

template<template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP>
struct _changeSign<GCD, CHKOP, false> {

    template<typename T>
    inline Rational<T, GCD, CHKOP> &operator() ( Rational<T, GCD, CHKOP> &r ) const {
        return r;
    };
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

}

}

namespace std {

/**
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
template<typename T, template<typename, bool, template<class, typename, bool> class> class GCD,
         template<class, typename, bool> class CHKOP> inline Commons::Math::Rational<T, GCD, CHKOP>
modf ( const Commons::Math::Rational<T, GCD, CHKOP> &__x,
       typename Commons::Math::Rational<T, GCD, CHKOP>::integer_type * __iptr ) {

    const typename Commons::Math::Rational<T, GCD, CHKOP>::mod_type &tmp ( __x.mod() );

    *__iptr = tmp.first;

    return tmp.second;
}

}

#endif /* COMMONS_MATH_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
