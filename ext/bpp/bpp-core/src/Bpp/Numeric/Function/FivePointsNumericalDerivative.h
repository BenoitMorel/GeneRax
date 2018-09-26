//
// File: FivePointsNumericalDerivative.h
// Created by: Julien Dutheil
// Created on: Thu Aug 17 15:00 2006
//

/*
   Copyright or © or Copr. CNRS, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for numerical calculus.

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

#ifndef _FIVEPOINTSNUMERICALDERIVATIVE_H_
#define _FIVEPOINTSNUMERICALDERIVATIVE_H_

#include "Functions.h"
#include "AbstractNumericalDerivative.h"

// From the STL:
#include <map>
#include <vector>
#include <string>

namespace bpp
{
/**
 * @brief Five points numerical derivative function wrapper.
 *
 * Numerical derivatives use three points to compute the derivatives.
 * @f$x_0@f$ is the focus point, @f$x_{-2} = x_0-2h@f$, @f$x_{-1} = x_0-h@f$, @f$x_{+1}=x_0+h@f$ and @f$x_{+2} = x_0+2h@f$
 * are other points, with function values @f$f_0@f$, @f$f_{-2}@f$, @f$f_{-1}@f$, @f$f_{+1}@f$ and @f$f_{+2}@f$ respectively.
 * The derivatives are then computed using the central formulas:
 * @f{eqnarray*}
 * \dfrac{\partial   f}{\partial x  } &=& \dfrac{f_{-2}-8f_{-1}+8f_{+1}-f_{+2}}{12h}\\
 * \dfrac{\partial^2 f}{\partial x^2} &=& \dfrac{-f_{-2}+16f_{-1}-30f_0+16f_{+1}-f_{+2}}{12h^2}\\
 * @f}
 * In case of border limit (when @f$x_{-2}@f$ or @f$x_{+2}@f$ are not computable),
 * the foreward and backward three points computations are performed, respectively:
 * @f{eqnarray*}
 * \dfrac{\partial   f}{\partial x  } &=& \dfrac{f_{+1}-f_0}{h}\\
 * \dfrac{\partial^2 f}{\partial x^2} &=& \dfrac{f_{+2}-2f_{+1}+f_0}{h^2}\\
 * @f}
 * and
 * @f{eqnarray*}
 * \dfrac{\partial   f}{\partial x  } &=& \dfrac{f_0-f_{-1}}{h}\\
 * \dfrac{\partial^2 f}{\partial x^2} &=& \dfrac{f_0-2f_{-1}+f_{-2}}{h^2}\\
 * @f}
 *
 * The @f$h@f$ parameter is computed in a parameter dependent manner:
 * @f$ h = x \times e@f$, with @f$x \neq 0@f$ being the current parameter value.
 * If @f$x = 0@f$, @f$h = e@f$.
 * Default value is provided for @f$e@f$ and corresponds to the _h field.
 * The default value works fine in most cases, but you may want to change it using the setInterval method.
 *
 * @warning cross second order derivatives are not implemented with the five points method.
 * @see AbstractNumericalDerivative
 */
class FivePointsNumericalDerivative :
  public AbstractNumericalDerivative
{
private:
  double f1_, f2_, f3_, f4_, f5_;

public:
  FivePointsNumericalDerivative(Function* function) :
    AbstractNumericalDerivative(function),
    f1_(),
    f2_(),
    f3_(),
    f4_(),
    f5_() {}
  FivePointsNumericalDerivative(DerivableFirstOrder* function) :
    AbstractNumericalDerivative(function),
    f1_(),
    f2_(),
    f3_(),
    f4_(),
    f5_() {}
  FivePointsNumericalDerivative(DerivableSecondOrder* function) :
    AbstractNumericalDerivative(function),
    f1_(),
    f2_(),
    f3_(),
    f4_(),
    f5_() {}
  virtual ~FivePointsNumericalDerivative() {}

  FivePointsNumericalDerivative* clone() const { return new FivePointsNumericalDerivative(*this); }

public:
  double getValue() const throw (Exception)
  {
    return f3_;
  }

  double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const
  throw (Exception)
  {
    throw Exception("Unimplemented cross derivative.");
  }

protected:
  void updateDerivatives(const ParameterList parameters)
  throw (ParameterNotFoundException, ConstraintException);
};
} // end of namespace bpp.

#endif // _FIVEPOINTSNUMERICALDERIVATIVE_H_

