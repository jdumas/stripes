#include <cmath>

#include "Vector.h"

namespace DDG
{
   template<typename URBG>
   Vector Vector :: randSphere( URBG && gen )
   {
      double z1 = unitRand(gen)*2.;
      double z2 = unitRand(gen);

      return Vector(
            std::sqrt( z1*(2.-z1)) * std::cos(2.*M_PI*z2),
            std::sqrt( z1*(2.-z1)) * std::sin(2.*M_PI*z2),
            1.-z1
            );
   }
}

