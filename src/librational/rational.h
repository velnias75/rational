/*
 * Copyright 2015 by Heiko Schäfer <heiko@rangun.de>
 *
 * This file is part of rational.
 *
 * NetMauMau is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation, either version 3 of
 * the License, or (at your option) any later version.
 *
 * NetMauMau is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with NetMauMau.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef COMMONS_MATH_RATIONAL_H
#define COMMONS_MATH_RATIONAL_H

#include <istream>
#include <ostream>
#include <limits>
#include <cmath>

namespace Commons {

namespace Math {

template<typename T>
class Rational {
    template<typename> friend class Rational;
public:
    Rational() : m_nom ( 0 ), m_denom ( 1 ) {}

    template<typename U>
    explicit Rational ( const Rational<U> &o ) : m_nom ( static_cast<T> ( o.m_nom ) ),
        m_denom ( static_cast<T> ( o.m_denom ) ) {}

    Rational ( T n, T d )  : m_nom ( n ), m_denom ( d ) {
        euclid ( *this );
    }

    template<typename FloatType>
    Rational ( FloatType f ) : m_nom ( 0 ), m_denom ( 1 ) {

        T p[2] = { 0, 1 };
        T q[2] = { 1, 0 };

        FloatType x = f;

        do {

            const T n = static_cast<T> ( std::floor ( x ) );
            x = static_cast<FloatType> ( 1 ) / ( x - static_cast<FloatType> ( n ) );

            m_nom = p[0] + n * p[1];
            p[0] = p[1];
            p[1] = m_nom;

            m_denom = q[0] + n * q[1];
            q[0] = q[1];
            q[1] = m_denom;

        } while ( ! ( std::abs ( static_cast<FloatType> ( m_nom ) /
                                 static_cast<FloatType> ( m_denom ) - f ) <
                      std::numeric_limits<FloatType>::epsilon() ) );
    }

    Rational &operator= ( const Rational& o ) {

        if ( this != &o ) {
            m_nom = o.m_nom;
            m_denom = o.m_denom;
        }

        return *this;
    }

    template<typename FloatType>
    inline Rational &operator= ( FloatType f ) {
        return ( *this = Rational ( f ) );
    }

    template<typename FloatType>
    inline operator FloatType() const {
        return static_cast<FloatType> ( m_nom ) /static_cast<FloatType> ( m_denom );
    }

    inline T nominator() const {
        return m_nom;
    }

    inline T denominator() const {
        return m_denom;
    }

    Rational &invert() {

        const T aux = m_denom;

        m_denom = m_nom;
        m_nom   = aux;

        return euclid ( *this );
    }

    inline Rational inv() const {
        return Rational ( *this ).invert();
    }

    Rational& operator+= ( const Rational& o ) {

        m_nom   = ( m_nom * o.m_denom ) + ( m_denom * o.m_nom );
        m_denom = m_denom * o.m_denom;

        return euclid ( *this );
    }

    template<typename FloatType>
    inline friend FloatType &operator+= ( FloatType &f, const Rational& o ) {
        return ( f = static_cast<FloatType> ( Rational ( f ) + o ) );
    }

    inline Rational operator+ ( const Rational& o ) const {
        return ( Rational ( *this ) += o );
    }

    template<typename FloatType>
    inline friend FloatType operator+ ( FloatType f, const Rational& o ) {
        return static_cast<FloatType> ( Rational ( f ) + o );
    }

    Rational& operator-= ( const Rational& o ) {

        m_nom   = ( m_nom * o.m_denom ) - ( m_denom * o.m_nom );
        m_denom = m_denom * o.m_denom;

        return euclid ( *this );
    }

    template<typename FloatType>
    inline friend FloatType &operator-= ( FloatType &f, const Rational& o ) {
        return ( f = static_cast<FloatType> ( Rational ( f ) - o ) );
    }

    inline Rational operator- ( const Rational& o ) const {
        return ( Rational ( *this ) -= o );
    }

    template<typename FloatType>
    inline friend FloatType operator- ( FloatType f, const Rational& o ) {
        return static_cast<FloatType> ( Rational ( f ) - o );
    }

    Rational& operator*= ( const Rational& o ) {

        m_nom *= o.m_nom;
        m_denom *= o.m_denom;

        return *this;
    }

    template<typename FloatType>
    inline friend FloatType &operator*= ( FloatType &f, const Rational& o ) {
        return ( f = static_cast<FloatType> ( Rational ( f ) * o ) );
    }

    inline Rational operator* ( const Rational& o ) const {
        return ( Rational ( *this ) *= o );
    }

    template<typename FloatType>
    inline friend FloatType operator* ( FloatType f, const Rational& o ) {
        return static_cast<FloatType> ( Rational ( f ) * o );
    }

    inline Rational& operator/= ( const Rational& o ) {
        return ( *this *= o.inv() );
    }

    template<typename FloatType>
    inline friend FloatType &operator/= ( FloatType &f, const Rational& o ) {
        return ( f = static_cast<FloatType> ( Rational ( f ) / o ) );
    }

    inline Rational operator/ ( const Rational& o ) const {
        return ( Rational ( *this ) /= o );
    }

    template<typename FloatType>
    inline friend FloatType operator/ ( FloatType f, const Rational& o ) {
        return static_cast<FloatType> ( Rational ( f ) / o );
    }

    inline bool operator== ( const Rational &o ) const {
        return ! ( ( *this < o ) || ( *this > o ) );
    }

    template<typename FloatType>
    inline friend bool operator== ( FloatType f, const Rational& o ) {
        return ( Rational ( f ) == o );
    }

    template<typename FloatType>
    inline friend bool operator== ( const Rational& o, FloatType f ) {
        return ( o == Rational ( f ) );
    }

    inline bool operator!= ( const Rational &o ) const {
        return ! ( *this == o );
    }

    template<typename FloatType>
    inline friend bool operator!= ( FloatType f, const Rational& o ) {
        return ! ( Rational ( f ) == o );
    }

    template<typename FloatType>
    inline friend bool operator!= ( const Rational& o, FloatType f ) {
        return ! ( o == Rational ( f ) );
    }

    inline bool operator< ( const Rational &o ) const {
        return ( m_denom * o.m_denom ) > 0 ? ( m_nom * o.m_denom ) < ( o.m_nom * m_denom )
               : ( o.m_nom * m_denom ) < ( m_nom * o.m_denom );
    }

    template<typename FloatType>
    inline friend bool operator< ( FloatType f, const Rational& o ) {
        return ( Rational ( f ) < o );
    }

    template<typename FloatType>
    inline friend bool operator< ( const Rational& o, FloatType f ) {
        return ( o < Rational ( f ) );
    }

    inline bool operator<= ( const Rational &o ) const {
        return ! ( o < *this );
    }

    template<typename FloatType>
    inline friend bool operator<= ( FloatType f, const Rational& o ) {
        return ! ( o < Rational ( f ) );
    }

    template<typename FloatType>
    inline friend bool operator<= ( const Rational& o, FloatType f ) {
        return ! ( Rational ( f ) < o );
    }

    inline bool operator> ( const Rational &o ) const {
        return o < *this;
    }

    template<typename FloatType>
    inline friend bool operator> ( FloatType f, const Rational& o ) {
        return o < Rational ( f );
    }

    template<typename FloatType>
    inline friend bool operator> ( const Rational& o, FloatType f ) {
        return Rational ( f ) < o;
    }

    inline bool operator>= ( const Rational &o ) const {
        return ! ( *this < o );
    }

    template<typename FloatType>
    inline friend bool operator>= ( FloatType f, const Rational& o ) {
        return ! ( Rational ( f ) < o );
    }

    template<typename FloatType>
    inline friend bool operator>= ( const Rational& o, FloatType f ) {
        return ! ( o < Rational ( f ) );
    }

private:
    Rational &euclid ( const Rational &o ) {

        T a = o.nominator(), b = o.denominator();

        while ( b ) {
            const T h = a % b;
            a = b;
            b = h;
        }

        m_nom /= a;
        m_denom /= a;

        return *this;
    }

private:
    T m_nom;
    T m_denom;
};

template<typename T>
std::ostream &operator<< ( std::ostream &o, const Rational<T> &r ) {
    return ( o << r.nominator() << "/" << r.denominator() );
}

template<typename T>
std::istream &operator>> ( std::istream &i, Rational<T> &r ) {

    double d;

    i >> d;
    r = d;

    return i;
}

}

}

#endif /* COMMONS_MATH_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 

