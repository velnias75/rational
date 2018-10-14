/*
 * Copyright 2018 by Heiko Schäfer <heiko@rangun.de>
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
 *
 * This header contains specialization for use of
 * unlogical as underlying storage type.
 *
 * @author Heiko Schäfer <heiko@rangun.de>
 * @copyright 2018 by Heiko Schäfer <heiko@rangun.de>
 *
 * @defgroup infint InfInt extensions
 *
 * The header `infint_rational.h` contains specializations especially for
 * [InfInt](http://sercantutar.github.io/infint/) as underlying
 * as underlying storage type.\n
 */

#ifndef COMMONS_MATH_UNLOGICAL_RATIONAL_H
#define COMMONS_MATH_UNLOGICAL_RATIONAL_H

#ifdef __EXCEPTIONS
/**
 * @ingroup infint
 * @def INFINT_RATIONAL_EXCEPTIONS
 *
 * Enables unlogical exceptions
 */
#ifndef UNLOGICAL_RATIONAL_EXCEPTIONS
#define UNLOGICAL_USE_EXCEPTIONS 1
#endif
#endif

#include <number.h>

#include "rational.h"

namespace std {
    
template<> struct numeric_limits<Commons::Lab::unlogical::Number<int64_t, 64> > {

    static const bool is_specialized = true;
    static const bool is_signed = true;
    static const bool is_integer = true;
    static const bool is_exact = false;

    static const Commons::Lab::unlogical::Number<int64_t, 64> min() {
        return Commons::Lab::unlogical::Number<int64_t, 64>(numeric_limits<int64_t>::min());
    }

    static const Commons::Lab::unlogical::Number<int64_t, 64> max() {
        return Commons::Lab::unlogical::Number<int64_t, 64>(numeric_limits<int64_t>::max());
    }
};

}

namespace Commons {

namespace Math {
    
#if 0

template<> struct TYPE_CONVERT<InfInt> {

    RATIONAL_NOCOPYASSIGN ( TYPE_CONVERT<InfInt> );

    explicit TYPE_CONVERT ( const InfInt &v ) : val ( v ) {}

    template<typename U> U convert() const;

private:
    const InfInt &val;
};

template<> inline InfInt TYPE_CONVERT<InfInt>::convert<InfInt>() const {
    return val;
}

template<> inline long double TYPE_CONVERT<InfInt>::convert<long double>() const {

    long double ld;

    ( std::istringstream ( val.toString() ) ) >> ld;

    return ld;
}

template<> inline double TYPE_CONVERT<InfInt>::convert<double>() const {

    double d;

    ( std::istringstream ( val.toString() ) ) >> d;

    return d;
}

template<> inline float TYPE_CONVERT<InfInt>::convert<float>() const {

    float f;

    ( std::istringstream ( val.toString() ) ) >> f;

    return f;
}

template<> inline InfInt TYPE_CONVERT<float>::convert<InfInt>() const {

    std::ostringstream os;

    os.precision ( std::numeric_limits<float>::digits );
    os << val;

    return os.str();
}

template<> inline InfInt TYPE_CONVERT<double>::convert<InfInt>() const {

    std::ostringstream os;

    os.precision ( std::numeric_limits<double>::digits );
    os << val;

    return os.str();
}

template<> inline InfInt TYPE_CONVERT<long double>::convert<InfInt>() const {

    std::ostringstream os;

    os.precision ( std::numeric_limits<long double>::digits );
    os << val;

    return os.str();
}

#endif /* #if 0 */

/**
 * @ingroup infint
 * @brief Rational class based on InfInt
 */
typedef Rational<Commons::Lab::unlogical::Number<int64_t, 64>, GCD_euclid> unlogical_rational;

template<> struct CFRationalTraits<Commons::Lab::unlogical::Number<int64_t, 64> > {
    typedef unlogical_rational rational_type;
};

}

}

#endif /* COMMONS_MATH_UNLOGICAL_RATIONAL_H */

// kate: indent-mode cstyle; indent-width 4; replace-tabs on; 
