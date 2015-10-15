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

#ifndef COMMONS_MATH_INFINT_RATIONAL_H
#define COMMONS_MATH_INFINT_RATIONAL_H

#ifdef COMMONS_MATH_GMP_RATIONAL_H
#error "infint_rational.h and gmp_rational.h cannot coexist in the same compilation unit"
#endif

#ifdef __EXCEPTIONS
#ifndef INFINT_RATIONAL_EXCEPTIONS
#define INFINT_USE_EXCEPTIONS 1
#endif
#endif

#include <InfInt.h>

#include "rational.h"

namespace std {

template<> struct numeric_limits<InfInt> {

     static const bool is_specialized = true;
     static const bool is_signed = true;
     static const bool is_integer = true;
     static const bool is_exact = false;

     static const InfInt min() {
          return InfInt();
     }

     static const InfInt max() {
          return InfInt();
     }
};

}

namespace Commons {

namespace Math {

template<> struct TYPE_CONVERT<long double> {

     inline explicit TYPE_CONVERT ( long double v ) : val ( v ) {}

     template<class U>
     inline U convert() const {

          std::ostringstream os;

          os.precision ( std::numeric_limits<long double>::digits );
          os << val;

          return os.str();
     }

private:
     long double val;
};

template<> inline long TYPE_CONVERT<long double>::convert<long>() const
{
     return std::floor ( val );
}

template<> struct TYPE_CONVERT<double> {

     inline explicit TYPE_CONVERT ( double v ) : val ( v ) {}

     template<class U>
     inline U convert() const {

          std::ostringstream os;

          os.precision ( std::numeric_limits<double>::digits );
          os << val;

          return os.str();
     }

private:
     double val;
};

template<> inline long TYPE_CONVERT<double>::convert<long>() const
{
     return std::floor ( val );
}

template<> struct TYPE_CONVERT<float> {

     inline explicit TYPE_CONVERT ( float v ) : val ( v ) {}

     template<class U>
     inline U convert() const {

          std::ostringstream os;

          os.precision ( std::numeric_limits<float>::digits );
          os << val;

          return os.str();
     }

private:
     float val;
};

template<> struct TYPE_CONVERT<InfInt> {

     inline explicit TYPE_CONVERT ( const InfInt &v ) : val ( v ) {}

     template<typename U> U convert() const;

private:
     const InfInt &val;
};

template<> inline InfInt TYPE_CONVERT<InfInt>::convert<InfInt>() const
{
     return val;
}

template<> inline long double TYPE_CONVERT<InfInt>::convert<long double>() const
{
     std::istringstream is ( val.toString() );

     long double ld;
     is >> ld;

     return ld;
}

template<> inline double TYPE_CONVERT<InfInt>::convert<double>() const
{
     std::istringstream is ( val.toString() );

     double d;
     is >> d;

     return d;
}

template<> inline float TYPE_CONVERT<InfInt>::convert<float>() const
{
     std::istringstream is ( val.toString() );

     float f;
     is >> f;

     return f;
}

typedef Rational<InfInt, GCD_euclid> infint_rational;

}

}

#endif /* COMMONS_MATH_INFINT_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 5; replace-tabs on; 
