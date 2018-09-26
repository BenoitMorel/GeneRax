//
// File: NumConstants.h
// Created by: Julien Dutheil
// Created on: Tue Feb 03 14:21 2009
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

  This software is a computer program whose purpose is to provide classes
  for numerical calculus. This file is part of the Bio++ project.

  This software is governed by the CeCILL  license under French law and
  abiding by the rules of distribution of free software.  You can  use, 
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info". 

  As a counterpart to the access to the source code and  rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty  and the software's author,  the holder of the
  economic rights,  and the successive licensors  have only  limited
  liability. 

  In this respect, the user's attention is drawn to the risks associated
  with loading,  using,  modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean  that it is complicated to manipulate,  and  that  also
  therefore means  that it is reserved for developers  and  experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or 
  data to be ensured and,  more generally, to use and operate it in the 
  same conditions as regards security. 

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _NUMCONSTANTS_H_
#define _NUMCONSTANTS_H_

#include <cmath>
#include <limits>

namespace bpp {

  /**
   * @brief this static class contains several useful constant values.
   *
   * This classe uses function in order to avoid the infamous "static initialization order fiasco".
   * C++0x solves this...
   */
  class NumConstants
  {

  public:
    /**
     * @name Golden ratio.
     *
     * The golden ratio, @f$\phi@f$ is equal to @f$\frac{1+\sqrt{5}}{2} = 1.6180339887498948482\ldots@f$.
     * We also define @f$R=\phi-1@f$ and @f$C = 1 - R@f$.
     * @{
     */
    static double GOLDEN_RATIO_PHI() { return (1. + sqrt(5.)) / 2.; }
    static double GOLDEN_RATIO_R() { return GOLDEN_RATIO_PHI() - 1.; }
    static double GOLDEN_RATIO_C() { return 1. - GOLDEN_RATIO_R(); }

    /** @} */

    static double MEGA() { return 1e6; }
    static double KILO() { return 1e3; }
    static double DECI() { return 1e-1; }
    static double CENTI() { return 1e-2; }
    static double MILLI() { return 1e-3; }
    static double MICRO() { return 1e-6; }
    static double NANO() { return 1e-9; }
    static double PICO() { return 1e-12; }

    static double SMALL() { return 1e-6; }
    static double TINY() { return 1e-12; }
    static double VERY_TINY() { return 1e-20; }
    static double VERY_BIG() { return static_cast<double>(1.7E+23); }
    
    /**
     * @name Define those constants in case they would not be available in stl/limits.
     *
     * @{
     */
    static double INF() { return std::numeric_limits<double>::has_infinity ? -log(0) : std::numeric_limits<double>::max(); }
    static double PINF() { return std::numeric_limits<double>::has_infinity ? -log(0) : std::numeric_limits<double>::max(); }
    static double MINF() { return std::numeric_limits<double>::has_infinity ? log(0) : std::numeric_limits<double>::min(); }
    static double NaN() { return NAN; }
    /** @} */

    static double PI() { return 3.141593; }
  };

}//end of namespace bpp.

#endif	//_NUMCONSTANTS_H_

