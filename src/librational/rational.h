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

template<template<typename, bool> class, bool>
struct _changeSign;

template<typename, template<typename, bool> class, bool>
struct _mod;

template<typename, template<typename, bool> class, bool>
struct _abs;

template<typename, bool>
struct GCD_stein;

template<typename, bool>
struct GCD_euclid;

template<typename, template<typename, bool> class, bool>
struct _lcm;

template<typename, template<typename, bool> class, typename, bool>
struct _approxFract;

template<typename T, template<typename, bool> class GCD = GCD_euclid>
class Rational {
    friend struct _changeSign<GCD, std::numeric_limits<T>::is_signed>;
    friend struct _mod<T, GCD, std::numeric_limits<T>::is_signed>;
    template<typename, bool> friend struct GCD_stein;
    template<typename, template<typename, bool> class, typename, bool> friend struct _approxFract;
public:
    typedef T integer_type;
    typedef typename _mod<integer_type, GCD,
            std::numeric_limits<integer_type>::is_signed>::pair_type mod_type;

    Rational() : m_numer (), m_denom ( 1 ) {}

    Rational ( const Rational &o ) : m_numer ( o.m_numer ), m_denom ( o.m_denom ) {}

    Rational ( const integer_type &n, const integer_type &d );

    Rational ( const integer_type &w, const integer_type &n, const integer_type &d ) : m_numer (),
        m_denom ( 1 ) {
        *this += Rational ( n, d ) += Rational ( w );
    }

    template<typename NumberType>
    Rational ( const NumberType &n );

    Rational &operator= ( const Rational& o ) {

        if ( this != &o ) {
            m_numer = o.m_numer;
            m_denom = o.m_denom;
        }

        return *this;
    }

    template<typename NumberType>
    inline Rational &operator= ( const NumberType &n ) {
        return ( *this = Rational ( n ) );
    }

    template<typename NumberType>
    inline operator NumberType() const {
        return static_cast<NumberType> ( m_numer ) / static_cast<NumberType> ( m_denom );
    }

    inline integer_type numerator() const throw() {
        return m_numer;
    }

    inline integer_type denominator() const throw() {
        return m_denom;
    }

    inline mod_type mod() const {
        return _mod<integer_type, GCD, std::numeric_limits<integer_type>::is_signed>() ( *this );
    }

    inline Rational abs() const {
        return _abs<integer_type, GCD, std::numeric_limits<integer_type>::is_signed>() ( *this );
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

    inline Rational operator+() const {
        return Rational ( *this );
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

    inline Rational operator-() const {
        Rational tmp ( *this );
        tmp.m_numer = -tmp.m_numer;
        return tmp;
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

    inline bool operator!() const {
        return m_numer == integer_type();
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
    Rational &gcm ( const Rational &o );

    static integer_type stein ( const integer_type &a, const integer_type &b ) {

        integer_type x ( a ), y ( b ), f = integer_type();

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

private:
    integer_type m_numer;
    integer_type m_denom;
};

template<typename T, template<typename, bool> class GCD>
Rational<T, GCD>::Rational ( const integer_type &n, const integer_type &d )  : m_numer ( n ),
    m_denom ( d ) {

    if ( m_denom == integer_type() ) throw std::runtime_error ( "denominator can't be null" );

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

    m_numer /= x;
    m_denom /= x;

    return _changeSign<GCD, std::numeric_limits<integer_type>::is_signed>() ( *this );
}

template<typename T, template<typename, bool> class GCD>
Rational<T, GCD>& Rational<T, GCD>::operator+= ( const Rational& o ) {

    if ( m_denom != o.m_denom ) {

        const integer_type &l ( _lcm<integer_type, GCD,
                                std::numeric_limits<integer_type>::is_signed>()
                                ( m_denom, o.m_denom ) );

        m_numer = ( ( l/m_denom ) * m_numer ) + ( ( l/o.m_denom ) * o.m_numer );
        m_denom = l;

    } else {
        m_numer = m_numer + o.m_numer;
    }

    return gcm ( *this );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline NumberType &operator+= ( NumberType &n, const Rational<T, GCD>& o ) {
    return ( n = Rational<T, GCD> ( n ) += o );
}

template<typename T, template<typename, bool> class GCD, typename NumberType>
inline Rational<T, GCD> operator+ ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o + Rational<T, GCD> ( n ) );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline Rational<T, GCD> operator+ ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) + o );
}

template<typename T, template<typename, bool> class GCD>
Rational<T, GCD>& Rational<T, GCD>::operator-= ( const Rational& o ) {

    if ( m_denom != o.m_denom ) {

        const integer_type &l ( _lcm<integer_type, GCD,
                                std::numeric_limits<integer_type>::is_signed>()
                                ( m_denom, o.m_denom ) );

        m_numer = ( ( l/m_denom ) * m_numer ) - ( ( l/o.m_denom ) * o.m_numer );
        m_denom = l;

    } else {
        m_numer = m_numer - o.m_numer;
    }

    return gcm ( *this );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline NumberType &operator-= ( NumberType &n, const Rational<T, GCD>& o ) {
    return ( n = Rational<T, GCD> ( n ) -= o );
}

template<typename T, template<typename, bool> class GCD, typename NumberType>
inline Rational<T, GCD> operator- ( const Rational<T, GCD>& o, const NumberType &n ) {
    return o - Rational<T, GCD> ( n );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline Rational<T, GCD> operator- ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) - o );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline NumberType &operator*= ( NumberType &n, const Rational<T, GCD>& o ) {
    return ( n = Rational<T, GCD> ( n ) *= o );
}

template<typename T, template<typename, bool> class GCD, typename NumberType>
inline Rational<T, GCD> operator* ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o * Rational<T, GCD> ( n ) );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline Rational<T, GCD> operator* ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) * o );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline NumberType &operator/= ( NumberType &n, const Rational<T, GCD>& o ) {
    return ( n = Rational<T, GCD> ( n ) /= o );
}

template<typename T, template<typename, bool> class GCD, typename NumberType>
inline Rational<T, GCD> operator/ ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o / Rational<T, GCD> ( n ) );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline Rational<T, GCD> operator/ ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) / o );
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

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline NumberType &operator%= ( NumberType &n, const Rational<T, GCD>& o ) {
    return ( n = Rational<T, GCD> ( n ) %= o );
}

template<typename T, template<typename, bool> class GCD, typename NumberType>
inline Rational<T, GCD> operator% ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o % Rational<T, GCD> ( n ) );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline Rational<T, GCD> operator% ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) % o );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator== ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) == o );
}

template<typename T, template<typename, bool> class GCD, typename NumberType>
inline bool operator== ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o == Rational<T, GCD> ( n ) );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator!= ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ! ( Rational<T, GCD> ( n ) == o );
}

template<typename T, template<typename, bool> class GCD, typename NumberType>
inline bool operator!= ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ! ( o == Rational<T, GCD> ( n ) );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator< ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ( Rational<T, GCD> ( n ) < o );
}

template<typename T, template<typename, bool> class GCD, typename NumberType>
inline bool operator< ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ( o < Rational<T, GCD> ( n ) );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator<= ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ! ( o < Rational<T, GCD> ( n ) );
}

template<typename T, template<typename, bool> class GCD, typename NumberType>
inline bool operator<= ( const Rational<T, GCD>& o, const NumberType &n ) {
    return ! ( Rational<T, GCD> ( n ) < o );
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator> ( const NumberType &n, const Rational<T, GCD>& o ) {
    return o < Rational<T, GCD> ( n );
}

template<typename T, template<typename, bool> class GCD, typename NumberType>
inline bool operator> ( const Rational<T, GCD>& o, const NumberType &n ) {
    return Rational<T, GCD> ( n ) < o;
}

template<typename NumberType, typename T, template<typename, bool> class GCD>
inline bool operator>= ( const NumberType &n, const Rational<T, GCD>& o ) {
    return ! ( Rational<T, GCD> ( n ) < o );
}

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

            const T n = static_cast<T> ( std::floor ( x ) );
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

template<typename T, bool>
struct GCD_euclid {

    inline T operator() ( const T &a, const T &b ) {

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
struct GCD_stein<T, true> {

    inline T operator() ( const T &a, const T &b ) {
        return Rational<T, GCD_stein::template GCD_stein>::stein ( a < T() ? -a : a,
                b < T() ? -b : b );
    }
};

template<typename T>
struct GCD_stein<T, false> {

    inline T operator() ( const T &a, const T &b ) {
        return Rational<T, GCD_stein::template GCD_stein>::stein ( a, b );
    }
};

template<typename T, template<typename, bool> class GCD>
struct _lcm<T, GCD, true> {

    inline T operator() ( const T &a, const T &b ) {

        const T &x ( a < T() ? -a : a ), &y ( b < T() ? -b : b );

        return ( static_cast<T> ( x ) / ( a ? GCD<T, false>() ( x, y ) : b ) ) *
               static_cast<T> ( y ) ;
    }
};

template<typename T, template<typename, bool> class GCD>
struct _lcm<T, GCD, false> {

    inline T operator() ( const T &a, const T &b ) {
        return ( a / ( a ? GCD<T, false>() ( a, b ) : b ) * b );
    }
};

template<template<typename, bool> class GCD>
struct _changeSign<GCD, true> {

    template<typename T>
    inline Rational<T, GCD> &operator() ( Rational<T, GCD> &r ) {

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
    inline Rational<T, GCD> &operator() ( Rational<T, GCD> &r ) {
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
