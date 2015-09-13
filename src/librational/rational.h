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

template<template<typename, bool> class, bool>
struct _changeSign;

template<typename, template<typename, bool> class, bool>
struct _mod;

template<typename, template<typename, bool> class, bool>
struct _abs;

/**
 * @brief Stein GCD algorithm implementation
 *
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 */
template<typename T, bool IsSigned>
struct GCD_stein;

/**
 * @brief Euclid GCD algorithm implementation
 *
 * @tparam T storage type
 * @tparam IsSigned specialization for @em signed or @em unsigned types
 */
template<typename T, bool IsSigned>
struct GCD_euclid;

template<typename, template<typename, bool> class, bool>
struct _lcm;

template<typename, template<typename, bool> class, typename, bool>
struct _approxFract;

/**
 * @brief %Rational (fraction) template class
 *
 * @note All %Rational objects are reduced
 *
 * @tparam T storage type
 * @tparam GCD GCD algorithm
 */
template<typename T, template<typename, bool> class GCD = GCD_euclid>
class Rational {
    friend struct _changeSign<GCD, std::numeric_limits<T>::is_signed>;
    friend struct _mod<T, GCD, std::numeric_limits<T>::is_signed>;
    template<typename, template<typename, bool> class, typename, bool> friend struct _approxFract;
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
    typedef typename _mod<integer_type, GCD,
            std::numeric_limits<integer_type>::is_signed>::pair_type mod_type;

    /**
     * @brief Creates a default (null) %Rational
     */
    Rational() : m_numer (), m_denom ( 1 ) {}

    /**
     * @brief Copy constructor
     *
     * @param[in] other the %Rational to copy
     */
    Rational ( const Rational &other ) : m_numer ( other.m_numer ), m_denom ( other.m_denom ) {}

    /**
     * @brief Creates a %Rational
     *
     * @param[in] numer the numerator
     * @param[in] denom the denominator
     */
    Rational ( const integer_type &numer, const integer_type &denom );

    /**
     * @brief Creates a inproper (mixed) %Rational
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
     * @brief Creates an approximated %Rational
     *
     * @tparam NumberType type of the number to approximate
     * @param[in] number the number to create an approximated %Rational of
     */
    template<typename NumberType>
    Rational ( const NumberType &number );

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
     * @brief assigns a Number
     *
     * The Number is approximated to a %Rational and then it gets assigned
     *
     * @tparam NumberType type of the number to approximate
     * @param[in] number the number to assign
     */
    template<typename NumberType>
    inline Rational &operator= ( const NumberType &number ) {
        return ( *this = Rational ( number ) );
    }

    /**
     * @brief convert to a Number
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
        return _mod<integer_type, GCD, std::numeric_limits<integer_type>::is_signed>() ( *this );
    }

    /**
     * @brief gets the absolute %Rational
     *
     * @return a copy of the absolute %Rational
     */
    inline Rational abs() const {
        return _abs<integer_type, GCD, std::numeric_limits<integer_type>::is_signed>() ( *this );
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
        return knuth_addSub ( other, std::plus<integer_type>() );
    }

    template<template<typename, bool> class U>
    inline friend Rational &operator+= ( Rational& x, const Rational<integer_type, U>& y ) {
        return ( x += Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief add a %Rational
     *
     * @param[in] other the %Rational to add
     *
     * @return a new %Rational
     */
    inline Rational operator+ ( const Rational& other ) const {
        return ( Rational ( *this ) += other );
    }

    template<template<typename, bool> class U>
    inline friend Rational operator+ ( const Rational& x, const Rational<integer_type, U>& y ) {
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
        m_numer += m_denom;
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
        return knuth_addSub ( other, std::minus<integer_type>() );
    }

    template<template<typename, bool> class U>
    inline friend Rational &operator-= ( Rational& x, const Rational<integer_type, U>& y ) {
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

    template<template<typename, bool> class U>
    inline friend Rational operator- ( const Rational& x, const Rational<integer_type, U>& y ) {
        return ( Rational ( x ) -= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief get a negated copy of the %Rational
     *
     * @return a negated copy of %Rational
     */
    inline Rational operator-() const {
        Rational tmp ( *this );
        tmp.m_numer = -tmp.m_numer;
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
        m_numer -= m_denom;
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

    template<template<typename, bool> class U>
    inline friend Rational &operator*= ( Rational& x, const Rational<integer_type, U>& y ) {
        return ( x *= Rational ( y.numerator(), y.denominator() ) );
    }

    /**
     * @brief multiply a %Rational
     *
     * @param[in] other the %Rational to multiply
     *
     * @return a new %Rational
     */
    inline Rational operator* ( const Rational& other ) const {
        return ( Rational ( *this ) *= other );
    }

    template<template<typename, bool> class U>
    inline friend Rational operator* ( const Rational& x, const Rational<integer_type, U>& y ) {
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

    template<template<typename, bool> class U>
    inline friend Rational &operator/= ( Rational& x, const Rational<integer_type, U>& y ) {
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

    template<template<typename, bool> class U>
    inline friend Rational operator/ ( const Rational& x, const Rational<integer_type, U>& y ) {
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

    template<template<typename, bool> class U>
    inline friend Rational &operator%= ( Rational& x, const Rational<integer_type, U>& y ) {
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

    template<template<typename, bool> class U>
    inline friend Rational operator% ( const Rational& x, const Rational<integer_type, U>& y ) {
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

    template<template<typename, bool> class U>
    inline friend bool operator== ( const Rational& x, const Rational<integer_type, U>& y ) {
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

    template<template<typename, bool> class U>
    inline friend bool operator!= ( const Rational& x, const Rational<integer_type, U>& y ) {
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

    template<template<typename, bool> class U>
    inline friend bool operator< ( const Rational& x, const Rational<integer_type, U>& y ) {
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

    template<template<typename, bool> class U>
    inline friend bool operator<= ( const Rational& x, const Rational<integer_type, U>& y ) {
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

    template<template<typename, bool> class U>
    inline friend bool operator> ( const Rational& x, const Rational<integer_type, U>& y ) {
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

    template<template<typename, bool> class U>
    inline friend bool operator>= ( const Rational& x, const Rational<integer_type, U>& y ) {
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
     * @param[in] mixed if @c true, than a mixed fraction is generated
     *
     * @return the string representation of %Rational
     */
    std::string str ( bool mixed = false ) const;

    friend std::ostream &operator<< ( std::ostream &o, const Rational &r ) {
        return ( o << r.str() );
    }

    friend std::istream &operator>> ( std::istream &i, Rational &r ) {

        double d;

        i >> d;
        r = d;

        return i;
    }

private:
    Rational &gcm ( const Rational &o );

    template<class Op>
    Rational &knuth_addSub ( const Rational &o, const Op &op );

private:
    integer_type m_numer;
    integer_type m_denom;
};

template<typename T, template<typename, bool> class GCD>
Rational<T, GCD>::Rational ( const integer_type &n, const integer_type &d )  : m_numer ( n ),
    m_denom ( d ) {

#ifdef __EXCEPTIONS
    if ( m_denom == integer_type() ) throw std::domain_error ( "denominator can't be null" );
#endif

    gcm ( *this );
}

template<typename T, template<typename, bool> class GCD> template<typename NumberType>
Rational<T, GCD>::Rational ( const NumberType &nt ) : m_numer ( static_cast<integer_type> ( nt ) ),
    m_denom ( 1 ) {

    _approxFract<integer_type, GCD, NumberType, ! ( std::numeric_limits<NumberType>::is_integer ||
            std::numeric_limits<NumberType>::is_exact ) >() ( *this, nt );
}

template<typename T, template<typename, bool> class GCD>
Rational<T, GCD> &Rational<T, GCD>::gcm ( const Rational &o ) {

    const integer_type &x ( o.m_numer ? GCD<T, std::numeric_limits<integer_type>::is_signed>()
                            ( o.m_numer, o.m_denom ) : o.m_denom );

    if ( x != static_cast<integer_type> ( 1 ) ) {
        m_numer /= x;
        m_denom /= x;
    }

    return _changeSign<GCD, std::numeric_limits<integer_type>::is_signed>() ( *this );
}

template<typename T, template<typename, bool> class GCD> template<class Op>
Rational<T, GCD> &Rational<T, GCD>::knuth_addSub ( const Rational<T, GCD> &o, const Op &op ) {

    const integer_type &d1 ( GCD<integer_type, std::numeric_limits<integer_type>::is_signed>()
                             ( m_denom, o.m_denom ) );

    if ( d1 == static_cast<integer_type> ( 1 ) ) {

        m_numer = op ( ( m_numer * o.m_denom ), ( m_denom * o.m_numer ) );
        m_denom *= o.m_denom;

    } else {

        const integer_type &t ( op ( m_numer * ( o.m_denom / d1 ),
                                     o.m_numer * ( m_denom / d1 ) ) );
        const integer_type &d2 ( GCD<integer_type, std::numeric_limits<integer_type>::is_signed>()
                                 ( t, d1 ) );
        m_numer = t / d2;
        m_denom = ( m_denom / d1 ) * ( o.m_denom / d2 );
    }

    return *this;
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline NumberType &operator+= ( NumberType &n, const Rational<T, GCD>& o ) {
    return ( n = Rational<T, GCD> ( n ) += o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 */
template<typename T, template<typename, bool> class GCD, typename NumberType>
inline Rational<T, GCD> operator+ ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o + Rational<T, GCD> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline Rational<T, GCD> operator+ ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) + o );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline NumberType &operator-= ( NumberType &n, const Rational<T, GCD>& o ) {
    return ( n = Rational<T, GCD> ( n ) -= o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 */
template<typename T, template<typename, bool> class GCD, typename NumberType>
inline Rational<T, GCD> operator- ( const Rational<T, GCD>& o, const NumberType &n ) {
    return o - Rational<T, GCD> ( n );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline Rational<T, GCD> operator- ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) - o );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline NumberType &operator*= ( NumberType &n, const Rational<T, GCD>& o ) {
    return ( n = Rational<T, GCD> ( n ) *= o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 */
template<typename T, template<typename, bool> class GCD, typename NumberType>
inline Rational<T, GCD> operator* ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o * Rational<T, GCD> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline Rational<T, GCD> operator* ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) * o );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline NumberType &operator/= ( NumberType &n, const Rational<T, GCD>& o ) {
    return ( n = Rational<T, GCD> ( n ) /= o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 */
template<typename T, template<typename, bool> class GCD, typename NumberType>
inline Rational<T, GCD> operator/ ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o / Rational<T, GCD> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline Rational<T, GCD> operator/ ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) / o );
}

template<typename T, template<typename, bool> class GCD>
Rational<T, GCD>& Rational<T, GCD>::operator*= ( const Rational& other ) {

    const integer_type &d1 ( GCD<integer_type,
                             std::numeric_limits<integer_type>::is_signed>() ( m_numer,
                                     other.m_denom ) );
    const integer_type &d2 ( GCD<integer_type,
                             std::numeric_limits<integer_type>::is_signed>() ( m_denom,
                                     other.m_numer ) );

    if ( ! ( d1 == static_cast<integer_type> ( 1 ) &&
             d2 == static_cast<integer_type> ( 1 ) ) ) {

        m_numer = ( m_numer / d1 ) * ( other.m_numer / d2 );
        m_denom = ( m_denom / d2 ) * ( other.m_denom / d1 );

    } else {
        m_numer *= other.m_numer;
        m_denom *= other.m_denom;
    }

    return *this;
}

template<typename T, template<typename, bool> class GCD>
Rational<T, GCD>& Rational<T, GCD>::operator%= ( const Rational& o ) {

    if ( m_denom != o.m_denom ) {

        const integer_type &l ( _lcm<integer_type, GCD,
                                std::numeric_limits<integer_type>::is_signed>()
                                ( m_denom, o.m_denom ) );

        const integer_type &a ( ( ( l/o.m_denom ) * o.m_numer ) );

        m_numer = ( ( ( l/m_denom ) * m_numer ) % a + a ) % a;
        m_denom = l;

    } else {
        m_numer = ( m_numer % o.m_numer + o.m_numer ) % o.m_numer;
    }

    return gcm ( *this );
}

template<typename T, template<typename, bool> class GCD>
std::string Rational<T, GCD>::str ( bool mixed ) const {

    std::ostringstream os;

    if ( mixed ) {

        const mod_type &p ( mod() );

        if ( p.first != integer_type() ) os << p.first << ' ';

        os << p.second.str ( false );

    } else {
        os << m_numer << '/' << m_denom;
    }

    return os.str();
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline NumberType &operator%= ( NumberType &n, const Rational<T, GCD>& o ) {
    return ( n = Rational<T, GCD> ( n ) %= o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 */
template<typename T, template<typename, bool> class GCD, typename NumberType>
inline Rational<T, GCD> operator% ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o % Rational<T, GCD> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline Rational<T, GCD> operator% ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) % o );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator== ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) == o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 */
template<typename T, template<typename, bool> class GCD, typename NumberType>
inline bool operator== ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o == Rational<T, GCD> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator!= ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ! ( Rational<T, GCD> ( n ) == o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 */
template<typename T, template<typename, bool> class GCD, typename NumberType>
inline bool operator!= ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ! ( o == Rational<T, GCD> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator< ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) < o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 */
template<typename T, template<typename, bool> class GCD, typename NumberType>
inline bool operator< ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o < Rational<T, GCD> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator<= ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ! ( o < Rational<T, GCD> ( n ) );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 */
template<typename T, template<typename, bool> class GCD, typename NumberType>
inline bool operator<= ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ! ( Rational<T, GCD> ( n ) < o );
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator> ( const NumberType &n, const Rational<T, GCD>& o ) {
    return o < Rational<T, GCD> ( n );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 */
template<typename T, template<typename, bool> class GCD, typename NumberType>
inline bool operator> ( const Rational<T, GCD>& o, const NumberType &n ) {
    return Rational<T, GCD> ( n ) < o;
}

/**
 * @relates Rational
 *
 * @tparam NumberType the number type
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 */
template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator>= ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ! ( Rational<T, GCD> ( n ) < o );
}

/**
 * @relates Rational
 *
 * @tparam T the storage type
 * @tparam GCD the GCD algorithm
 * @tparam NumberType the number type
 */
template<typename T, template<typename, bool> class GCD, typename NumberType>
inline bool operator>= ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ! ( o < Rational<T, GCD> ( n ) );
}

template<typename T, template<typename, bool> class GCD, typename NumberType>
struct _approxFract<T, GCD, NumberType, true> {

    inline void operator() ( Rational<T, GCD> &r, const NumberType &nt ) const {

        T p[2] = { T(), 1 };
        T q[2] = { 1, T() };

        NumberType x ( nt );

        while ( ! ( abs ( static_cast<NumberType> ( r.m_numer ) /
                          static_cast<NumberType> ( r.m_denom ) - nt ) <
                    std::numeric_limits<NumberType>::epsilon() ) ) {

            const T &n ( static_cast<T> ( std::floor ( x ) ) );
            x = static_cast<NumberType> ( 1 ) / ( x - static_cast<NumberType> ( n ) );

            r.m_numer = p[0] + n * p[1];
            p[0] = p[1];
            p[1] = r.m_numer;

            r.m_denom = q[0] + n * q[1];
            q[0] = q[1];
            q[1] = r.m_denom;
        }
    }

private:
    inline NumberType abs ( const NumberType &nt ) const {
        return nt < NumberType() ? -nt : nt;
    }
};

template<typename T, template<typename, bool> class GCD, typename NumberType>
struct _approxFract<T, GCD, NumberType, false> {
    inline void operator() ( const Rational<T, GCD> &, const NumberType & ) const {}
};

template<typename T, template<typename, bool> class GCD>
struct _abs<T, GCD, true> {

    inline Rational<T, GCD> operator() ( const Rational<T, GCD> &r ) const {
        return r.numerator() < T() ? -r : r;
    }
};

template<typename T, template<typename, bool> class GCD>
struct _abs<T, GCD, false> {

    inline Rational<T, GCD> operator() ( const Rational<T, GCD> &r ) const {
        return r;
    }
};

template<typename T, template<typename, bool> class GCD>
struct _mod<T, GCD, true> {

    typedef std::pair<T, Rational<T, GCD> > pair_type;

    inline pair_type operator() ( const Rational<T, GCD> &r ) const {

        const Rational<T, GCD> &h ( Rational<T, GCD> ( ( r.m_numer % r.m_denom ), r.m_denom ) );

        return std::make_pair ( r.m_numer/r.m_denom, r.m_numer < T() ? -h : h );
    }
};

template<typename T, template<typename, bool> class GCD>
struct _mod<T, GCD, false> {

    typedef std::pair<T, Rational<T, GCD> > pair_type;

    inline pair_type operator() ( const Rational<T, GCD> &r ) const {
        return std::make_pair ( r.m_numer/r.m_denom, Rational<T, GCD> ( ( r.m_numer % r.m_denom ),
                                r.m_denom ) );
    }
};

template<typename T>
struct GCD_euclid<T, false> {

    inline T operator() ( const T &a, const T &b ) const {

        T x ( a ), y ( b );

        // while ( y ) { const integer_type &h ( x % y ); x = y; y = h; }

        while ( y ) {

            x %= y;
            y ^= x;
            x ^= y;
            y ^= x;
        }

        return x;
    }
};

template<typename T>
struct GCD_euclid<T, true> {

    inline T operator() ( const T &a, const T &b ) const {
        const T &h ( GCD_euclid<T, false>() ( a, b ) );
        return h < T() ? -h : h;
    }
};

template<typename T>
struct GCD_stein<T, false> {

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

template<typename T>
struct GCD_stein<T, true> {

    inline T operator() ( const T &a, const T &b ) const {
        return GCD_stein<T, false>() ( a < T() ? -a : a, b < T() ? -b : b );
    }
};

template<typename T, template<typename, bool> class GCD>
struct _lcm<T, GCD, true> {

    inline T operator() ( const T &a, const T &b ) const {

        const T &x ( a < T() ? -a : a ), &y ( b < T() ? -b : b );

        return ( static_cast<T> ( x ) / ( a ? GCD<T, false>() ( x, y ) : b ) ) *
               static_cast<T> ( y ) ;
    }
};

template<typename T, template<typename, bool> class GCD>
struct _lcm<T, GCD, false> {

    inline T operator() ( const T &a, const T &b ) const {
        return ( a / ( a ? GCD<T, false>() ( a, b ) : b ) * b );
    }
};

template<template<typename, bool> class GCD>
struct _changeSign<GCD, true> {

    template<typename T>
    inline Rational<T, GCD> &operator() ( Rational<T, GCD> &r ) const {

        if ( r.m_denom < T() ) {
            r.m_numer = -r.m_numer;
            r.m_denom = -r.m_denom;
        }

        return r;
    };
};

template<template<typename, bool> class GCD>
struct _changeSign<GCD, false> {

    template<typename T>
    inline Rational<T, GCD> &operator() ( Rational<T, GCD> &r ) const {
        return r;
    };
};

}

}

namespace std {

template<typename T, template<typename, bool> class GCD>
inline Commons::Math::Rational<T, GCD> modf ( const Commons::Math::Rational<T, GCD> &__x,
        typename Commons::Math::Rational<T, GCD>::integer_type * __iptr ) {

    const typename Commons::Math::Rational<T, GCD>::mod_type &tmp ( __x.mod() );

    *__iptr = tmp.first;

    return tmp.second;
}

}

#endif /* COMMONS_MATH_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
