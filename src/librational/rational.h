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

#ifndef COMMONS_MATH_RATIONAL_H
#define COMMONS_MATH_RATIONAL_H

#include <stdexcept>
#include <sstream>
#include <limits>
#include <cmath>

namespace Commons {

namespace Math {

template<bool>
struct _changeSign;

template<typename, bool>
struct _mod;

template<typename T>
class Rational {
    template<typename> friend class Rational;
    template<bool> friend struct _changeSign;
    template<typename, bool> friend struct _mod;
public:
    typedef T integer_type;
    typedef typename _mod<integer_type,
            std::numeric_limits<integer_type>::is_signed>::pair_type mod_type;

    Rational() : m_numer ( 0 ), m_denom ( 1 ) {}

    Rational ( const Rational &o ) : m_numer ( o.m_numer ), m_denom ( o.m_denom ) {}

    Rational ( const integer_type &n, const integer_type &d );

    template<typename FloatType>
    Rational ( const FloatType &f );

    Rational &operator= ( const Rational& o ) {

        if ( this != &o ) {
            m_numer = o.m_numer;
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
        return static_cast<FloatType> ( m_numer ) / static_cast<FloatType> ( m_denom );
    }

    inline integer_type numerator() const throw() {
        return m_numer;
    }

    inline integer_type denominator() const throw() {
        return m_denom;
    }

    inline mod_type mod() const {
        return _mod<integer_type, std::numeric_limits<integer_type>::is_signed>() ( *this );
    }

    Rational &invert() {

        using namespace std;
        swap ( m_numer, m_denom );

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
        m_numer += m_denom;
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
        m_numer -= m_denom;
        return gcm ( *this );
    }

    inline Rational operator-- ( int ) {
        Rational tmp ( *this );
        --*this;
        return tmp;
    }

    Rational& operator*= ( const Rational& o ) {

        m_numer *= o.m_numer;
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
        return /* ( m_denom * o.m_denom ) > 0 ? */ ( m_numer * o.m_denom ) < ( o.m_numer * m_denom )
                /* : ( o.m_numer * m_denom ) < ( m_numer * o.m_denom ) */;
        // denom can NEVER be zero!
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
    integer_type m_numer;
    integer_type m_denom;
};

template<typename T>
Rational<T>::Rational ( const integer_type &n, const integer_type &d )  : m_numer ( n ),
    m_denom ( d ) {

    if ( m_denom == integer_type() ) throw std::runtime_error ( "denominator can't be null" );

    gcm ( *this );
}

template<typename T> template<typename FloatType>
Rational<T>::Rational ( const FloatType &f ) : m_numer ( static_cast<integer_type> ( f ) ),
    m_denom ( 1 ) {

    if ( ! ( std::numeric_limits<FloatType>::is_integer ||
             std::numeric_limits<FloatType>::is_exact ) ) {

        integer_type p[2] = { 0, 1 };
        integer_type q[2] = { 1, 0 };

        FloatType x ( f );

        while ( ! ( std::abs ( static_cast<FloatType> ( m_numer ) /
                               static_cast<FloatType> ( m_denom ) - f ) <
                    std::numeric_limits<FloatType>::epsilon() ) ) {

            const integer_type n = static_cast<integer_type> ( std::floor ( x ) );
            x = static_cast<FloatType> ( 1 ) / ( x - static_cast<FloatType> ( n ) );

            m_numer = p[0] + n * p[1];
            p[0] = p[1];
            p[1] = m_numer;

            m_denom = q[0] + n * q[1];
            q[0] = q[1];
            q[1] = m_denom;

        }
    }
}

template<typename T>
Rational<T> &Rational<T>::gcm ( const Rational &o ) {

    const integer_type &x ( o.m_numer ? euclid ( o.m_numer, o.m_denom ) : o.m_denom );

    m_numer /= x;
    m_denom /= x;

    return _changeSign<std::numeric_limits<integer_type>::is_signed>() ( *this );
}

template<typename T>
Rational<T>& Rational<T>::operator+= ( const Rational& o ) {

    if ( m_denom != o.m_denom ) {

        const integer_type &l ( lcm ( m_denom, o.m_denom ) );

        m_numer = ( ( l/m_denom ) * m_numer ) + ( ( l/o.m_denom ) * o.m_numer );
        m_denom = l;

    } else {
        m_numer = m_numer + o.m_numer;
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

        m_numer = ( ( l/m_denom ) * m_numer ) - ( ( l/o.m_denom ) * o.m_numer );
        m_denom = l;

    } else {
        m_numer = m_numer - o.m_numer;
    }

    return gcm ( *this );
}

template<typename T, typename FloatType>
inline FloatType &operator-= ( FloatType &f, const Rational<T>& o ) {
    return ( f = Rational<T> ( f ) -= o );
}

template<typename T, typename FloatType>
inline Rational<T> operator- ( const Rational<T>& o, const FloatType &f ) {
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
inline Rational<T> operator/ ( const Rational<T>& o, const FloatType &f ) {
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
        const integer_type &a ( ( ( l/o.m_denom ) * o.m_numer ) );

        m_numer = ( ( ( l/m_denom ) * m_numer ) % a + a ) % a;
        m_denom = l;

    } else {
        m_numer = ( m_numer % o.m_numer + o.m_numer ) % o.m_numer;
    }

    return gcm ( *this );
}

template<typename T>
std::string Rational<T>::str ( bool mixed ) const {

    std::ostringstream os;

    if ( mixed ) {

        const std::pair<integer_type, Rational<T> > &p ( mod() );

        if ( p.first != 0 ) os << p.first << ' ';

        os  << p.second.str ( false );

    } else {
        os << m_numer << '/' << m_denom;
    }

    return os.str();
}

template<typename T, typename FloatType>
inline FloatType &operator%= ( FloatType &f, const Rational<T>& o ) {
    return ( f = Rational<T> ( f ) %= o );
}

template<typename T, typename FloatType>
inline Rational<T> operator% ( const Rational<T>& o, const FloatType &f ) {
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

template<typename T>
struct _mod<T, true> {

    typedef std::pair<T, Rational<T> > pair_type;

    inline pair_type operator() ( const Rational<T> &r ) const {
        return std::make_pair ( r.m_numer/r.m_denom, Rational<T> ( ( r.m_numer % r.m_denom ) *
                                ( r.m_numer < T() ? static_cast<T> ( -1 ) : static_cast<T> ( 1 ) ),
                                r.m_denom ) );
    }
};

template<typename T>
struct _mod<T, false> {

    typedef std::pair<T, Rational<T> > pair_type;

    inline pair_type operator() ( const Rational<T> &r ) const {
        return std::make_pair ( r.m_numer/r.m_denom, Rational<T> ( ( r.m_numer % r.m_denom ),
                                r.m_denom ) );
    }
};

template<>
struct _changeSign<true> {

    template<typename T>
    inline Rational<T> &operator() ( Rational<T> &r ) {

        if ( r.m_denom < T() ) {
            r.m_numer *= static_cast<T> ( -1 );
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
