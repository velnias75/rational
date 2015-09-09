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

#include <stdexcept>
#include <istream>
#include <ostream>
#include <limits>
#include <cmath>

namespace Commons {

namespace Math {

template<bool>
struct _changeSign;

template<typename T>
class Rational {
    template<typename> friend class Rational;
    template<bool> friend struct _changeSign;
public:
    typedef T integer_type;

    Rational() : m_nom ( 0 ), m_denom ( 1 ) {}

    Rational ( const Rational &o ) : m_nom ( o.m_nom ), m_denom ( o.m_denom ) {}

    Rational ( const integer_type &n, const integer_type &d )  : m_nom ( n ), m_denom ( d ) {

        if ( m_denom == integer_type() ) throw std::runtime_error ( "denominator can't be null" );

        gcm ( *this );
    }

    template<typename FloatType>
    Rational ( const FloatType &f );

    Rational &operator= ( const Rational& o ) {

        if ( this != &o ) {
            m_nom = o.m_nom;
            m_denom = o.m_denom;
        }

        return *this;
    }

    template<typename FloatType>
    inline Rational &operator= ( const FloatType &f ) {
        return ( *this = Rational ( f ) );
    }

    template<typename FloatType>
    inline operator FloatType() const {
        return static_cast<FloatType> ( m_nom ) / static_cast<FloatType> ( m_denom );
    }

    inline integer_type nominator() const {
        return m_nom;
    }

    inline integer_type denominator() const {
        return m_denom;
    }

    Rational &invert() {

        using namespace std;
        swap ( m_nom, m_denom );

        if ( m_denom == integer_type() ) throw std::runtime_error ( "division by zero" );

        return *this;
    }

    inline Rational inverse() const {
        return Rational ( *this ).invert();
    }

    Rational& operator+= ( const Rational& o );

    inline Rational operator+ ( const Rational& o ) const {
        return ( Rational ( *this ) += o );
    }

    inline Rational& operator++() {
        m_nom += m_denom;
        return gcm ( *this );
    }

    inline Rational operator++ ( int ) {
        Rational tmp ( *this );
        ++*this;
        return tmp;
    }

    Rational& operator-= ( const Rational& o );

    inline Rational operator- ( const Rational& o ) const {
        return ( Rational ( *this ) -= o );
    }

    inline Rational& operator--() {
        m_nom -= m_denom;
        return gcm ( *this );
    }

    inline Rational operator-- ( int ) {
        Rational tmp ( *this );
        --*this;
        return tmp;
    }

    Rational& operator*= ( const Rational& o ) {

        m_nom *= o.m_nom;
        m_denom *= o.m_denom;

        return *this;
    }

    inline Rational operator* ( const Rational& o ) const {
        return ( Rational ( *this ) *= o );
    }

    inline Rational& operator/= ( const Rational& o ) {
        return ( *this *= o.inverse() );
    }

    inline Rational operator/ ( const Rational& o ) const {
        return ( Rational ( *this ) /= o );
    }

    Rational& operator%= ( const Rational& o );

    inline Rational operator% ( const Rational& o ) const {
        return ( Rational ( *this ) %= o );
    }

    inline bool operator== ( const Rational &o ) const {
        return ! ( ( *this < o ) || ( *this > o ) );
    }

    inline bool operator!= ( const Rational &o ) const {
        return ! ( *this == o );
    }

    inline bool operator< ( const Rational &o ) const {
        return /* ( m_denom * o.m_denom ) > 0 ? */ ( m_nom * o.m_denom ) < ( o.m_nom * m_denom )
               /* : ( o.m_nom * m_denom ) < ( m_nom * o.m_denom ) */; // denom can NEVER be zero!
    }

    inline bool operator<= ( const Rational &o ) const {
        return ! ( o < *this );
    }

    inline bool operator> ( const Rational &o ) const {
        return o < *this;
    }

    inline bool operator>= ( const Rational &o ) const {
        return ! ( *this < o );
    }

    friend std::ostream &operator<< ( std::ostream &o, const Rational &r ) {
        return ( o << r.m_nom << "/" << r.m_denom );
    }

    friend std::istream &operator>> ( std::istream &i, Rational &r ) {

        double d;

        i >> d;
        r = d;

        return i;
    }

private:
    inline static integer_type lcm ( const integer_type &a, const integer_type &b ) {
        return std::numeric_limits<integer_type>::is_signed ?
               ( static_cast<integer_type> ( std::abs ( a ) ) / ( a ? euclid ( a, b ) : b ) ) *
               static_cast<integer_type> ( std::abs ( b ) ) : a / ( a ? euclid ( a, b ) : b ) * b ;
    }

    Rational &gcm ( const Rational &o );

    static integer_type euclid ( const integer_type &a, const integer_type &b ) {

        integer_type x ( a ), y ( b );

        // while ( y ) { const integer_type &h ( x % y ); x = y; y = h; }

        while ( y ) {
            x %= y;
            y ^= x;
            x ^= y;
            y ^= x;
        }

        return x;
    }

private:
    integer_type m_nom;
    integer_type m_denom;
};

template<typename T> template<typename FloatType>
Rational<T>::Rational ( const FloatType &f ) : m_nom ( 0 ), m_denom ( 1 ) {

    if ( !std::numeric_limits<FloatType>::is_exact ) {

        integer_type p[2] = { 0, 1 };
        integer_type q[2] = { 1, 0 };

        FloatType x ( f );

        do {

            const integer_type n = static_cast<integer_type> ( std::floor ( x ) );
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
    } else {
        m_nom = static_cast<integer_type> ( f );
        m_denom = static_cast<integer_type> ( 1 );
    }
}

template<typename T>
Rational<T> &Rational<T>::gcm ( const Rational &o ) {

    const integer_type &x ( o.m_nom ? euclid ( o.m_nom, o.m_denom ) : o.m_denom );

    m_nom /= x;
    m_denom /= x;

    return _changeSign<std::numeric_limits<integer_type>::is_signed>() ( *this );
}

template<typename T>
Rational<T>& Rational<T>::operator+= ( const Rational& o ) {

    if ( m_denom != o.m_denom ) {

        const integer_type &l ( lcm ( m_denom, o.m_denom ) );

        m_nom = ( ( l/m_denom ) * m_nom ) + ( ( l/o.m_denom ) * o.m_nom );
        m_denom = l;

    } else {
        m_nom = m_nom + o.m_nom;
    }

    return gcm ( *this );
}

template<typename T, typename FloatType>
inline FloatType &operator+= ( FloatType &f, const Rational<T>& o ) {
    return ( f = Rational<T> ( f ) += o );
}

template<typename T, typename FloatType>
inline Rational<T> operator+ ( const Rational<T>& o, const FloatType &f ) {
    return ( o + Rational<T> ( f ) );
}

template<typename T, typename FloatType>
inline Rational<T> operator+ ( const FloatType &f, const Rational<T>& o ) {
    return ( Rational<T> ( f ) + o );
}

template<typename T>
Rational<T>& Rational<T>::operator-= ( const Rational& o ) {

    if ( m_denom != o.m_denom ) {

        const integer_type &l ( lcm ( m_denom, o.m_denom ) );

        m_nom = ( ( l/m_denom ) * m_nom ) - ( ( l/o.m_denom ) * o.m_nom );
        m_denom = l;

    } else {
        m_nom = m_nom - o.m_nom;
    }

    return gcm ( *this );
}

template<typename T, typename FloatType>
inline FloatType &operator-= ( FloatType &f, const Rational<T>& o ) {
    return ( f = Rational<T> ( f ) -= o );
}

template<typename T, typename FloatType>
inline FloatType operator- ( const Rational<T>& o, const FloatType &f ) {
    return o - Rational<T> ( f );
}

template<typename T, typename FloatType>
inline Rational<T> operator- ( const FloatType &f, const Rational<T>& o ) {
    return ( Rational<T> ( f ) - o );
}

template<typename T, typename FloatType>
inline FloatType &operator*= ( FloatType &f, const Rational<T>& o ) {
    return ( f = Rational<T> ( f ) *= o );
}

template<typename T, typename FloatType>
inline Rational<T> operator* ( const Rational<T>& o, const FloatType &f ) {
    return ( o * Rational<T> ( f ) );
}

template<typename T, typename FloatType>
inline Rational<T> operator* ( const FloatType &f, const Rational<T>& o ) {
    return ( Rational<T> ( f ) * o );
}

template<typename T, typename FloatType>
inline FloatType &operator/= ( FloatType &f, const Rational<T>& o ) {
    return ( f = Rational<T> ( f ) /= o );
}

template<typename T, typename FloatType>
inline FloatType operator/ ( const Rational<T>& o, const FloatType &f ) {
    return ( o / Rational<T> ( f ) );
}

template<typename T, typename FloatType>
inline Rational<T> operator/ ( const FloatType &f, const Rational<T>& o ) {
    return ( Rational<T> ( f ) / o );
}

template<typename T>
Rational<T>& Rational<T>::operator%= ( const Rational& o ) {

    if ( m_denom != o.m_denom ) {

        const integer_type &l ( lcm ( m_denom, o.m_denom ) );
        const integer_type &a ( ( ( l/o.m_denom ) * o.m_nom ) );

        m_nom = ( ( ( l/m_denom ) * m_nom ) % a + a ) % a;
        m_denom = l;

    } else {
        m_nom = ( m_nom % o.m_nom + o.m_nom ) % o.m_nom;
    }

    return gcm ( *this );
}

template<typename T, typename FloatType>
inline FloatType &operator%= ( FloatType &f, const Rational<T>& o ) {
    return ( f = Rational<T> ( f ) %= o );
}

template<typename T, typename FloatType>
inline FloatType operator% ( const Rational<T>& o, const FloatType &f ) {
    return ( o % Rational<T> ( f ) );
}

template<typename T, typename FloatType>
inline Rational<T> operator% ( const FloatType &f, const Rational<T>& o ) {
    return ( Rational<T> ( f ) % o );
}

template<typename T, typename FloatType>
inline bool operator== ( const FloatType &f, const Rational<T>& o ) {
    return ( Rational<T> ( f ) == o );
}

template<typename T, typename FloatType>
inline bool operator== ( const Rational<T>& o, const FloatType &f ) {
    return ( o == Rational<T> ( f ) );
}

template<typename T, typename FloatType>
inline bool operator!= ( const FloatType &f, const Rational<T>& o ) {
    return ! ( Rational<T> ( f ) == o );
}

template<typename T, typename FloatType>
inline bool operator!= ( const Rational<T>& o, const FloatType &f ) {
    return ! ( o == Rational<T> ( f ) );
}

template<typename T, typename FloatType>
inline bool operator< ( const FloatType &f, const Rational<T>& o ) {
    return ( Rational<T> ( f ) < o );
}

template<typename T, typename FloatType>
inline bool operator< ( const Rational<T>& o, const FloatType &f ) {
    return ( o < Rational<T> ( f ) );
}

template<typename T, typename FloatType>
inline bool operator<= ( const FloatType &f, const Rational<T>& o ) {
    return ! ( o < Rational<T> ( f ) );
}

template<typename T, typename FloatType>
inline bool operator<= ( const Rational<T>& o, const FloatType &f ) {
    return ! ( Rational<T> ( f ) < o );
}

template<typename T, typename FloatType>
inline bool operator> ( const FloatType &f, const Rational<T>& o ) {
    return o < Rational<T> ( f );
}

template<typename T, typename FloatType>
inline bool operator> ( const Rational<T>& o, const FloatType &f ) {
    return Rational<T> ( f ) < o;
}

template<typename T, typename FloatType>
inline bool operator>= ( const FloatType &f, const Rational<T>& o ) {
    return ! ( Rational<T> ( f ) < o );
}

template<typename T, typename FloatType>
inline bool operator>= ( const Rational<T>& o, const FloatType &f ) {
    return ! ( o < Rational<T> ( f ) );
}

template<>
struct _changeSign<true> {

    template<typename T>
    inline Rational<T> &operator() ( Rational<T> &r ) {

        if ( r.m_denom < 0 ) {
            r.m_nom *= static_cast<T> ( -1 );
            r.m_denom *= static_cast<T> ( -1 );
        }

        return r;
    };
};

template<>
struct _changeSign<false> {

    template<typename T>
    inline Rational<T> &operator() ( Rational<T> &r ) {
        return r;
    };
};

}

}

#endif /* COMMONS_MATH_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on;
