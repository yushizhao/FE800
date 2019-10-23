/* -*- mode: c++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */

/*
 Copyright (C) 2003 Ferdinando Ametrano
 Copyright (C) 2010 Kakhkhor Abdijalilov

 This file is part of QuantLib, a free-software/open-source library
 for financial quantitative analysts and developers - http://quantlib.org/

 QuantLib is free software: you can redistribute it and/or modify it
 under the terms of the QuantLib license.  You should have received a
 copy of the license along with this program; if not, please email
 <quantlib-dev@lists.sf.net>. The license is also available online at
 <http://quantlib.org/license.shtml>.

 This program is distributed in the hope that it will be useful, but WITHOUT
 ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 FOR A PARTICULAR PURPOSE.  See the license for more details.
*/

/*! \file mt19937uniformrng.hpp
    \brief Mersenne Twister uniform random number generator
*/
#include <ql/quantlib.hpp>
#include <fstream>
#include <stdint.h> 


using namespace QuantLib;

class TRNG {     
      
      public:
        static std::ifstream binary;
        uint32_t a;
        int size;
        typedef Sample<Real> sample_type;
        /*! if the given seed is 0, a random seed will be chosen
            based on clock() */
        explicit TRNG(uint32_t seed = 0);

        sample_type next() const {
          binary.read((char*)&a, size);
          return sample_type(a / 65536.0 / 65536.0, 1.0); 
        };

        // explicit TRNG(
        //                              const std::vector<unsigned long>& seeds);
        /*! returns a sample with weight 1.0 containing a random number
            in the (0.0, 1.0) interval  */
        // sample_type next() const { return sample_type(nextReal(),1.0); }
        //! return a random number in the (0.0, 1.0)-interval
        // Real nextReal() const {
        //     return (Real(nextInt32()) + 0.5)/4294967296.0;
        // }
        //! return a random integer in the [0,0xffffffff]-interval
        // unsigned long nextInt32() const  {
        //     if (mti==N)
        //         twist(); /* generate N words at a time */

        //     unsigned long y = mt[mti++];

        //      Tempering 
        //     y ^= (y >> 11);
        //     y ^= (y << 7) & 0x9d2c5680UL;
        //     y ^= (y << 15) & 0xefc60000UL;
        //     y ^= (y >> 18);
        //     return y;
        // }
    };

std::ifstream TRNG::binary("ANU");

TRNG::TRNG(uint32_t seed) {
    size = sizeof(seed);
}
