//
// File: TwoPointsNumericalDerivative.h
// Created by: Julien Dutheil
// Created on: Mon May 28 10:33 2007
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

#ifndef _TWOPOINTSNUMERICALDERIVATIVE_H_
#define _TWOPOINTSNUMERICALDERIVATIVE_H_

#include "Functions.h"
#include "AbstractNumericalDerivative.h"

// From the STL:
#include <map>
#include <vector>
#include <string>

namespace bpp
{
/**
 * @brief Two points numerical derivative function wrapper.
 *
 * Numerical derivatives use two points to compute the (first order) derivatives.
 * @f$x_0@f$ is the focus point and @f$x_{+1} = x_0+h@f$ (or @f$x_0-h@f$, if constrained on the right).
 * Corresponding function values are @f$f_0@f$, @f$f_{-1}@f$ and @f$f_{+1/-1}@f$ respectively.
 * The derivatives are then computed using the central formulas:
 * @f{eqnarray*}
 * \dfrac{\partial   f}{\partial x  } &=& \dfrac{f_{+1}-f_{0}}{h}\mathrm{\ or}\\
 * \dfrac{\partial   f}{\partial x  } &=& \dfrac{f_{0}-f_{-1}}{h}\\
 * @f}
 * This class does not allow computation of second order derivatives.
 *
 * The @f$h@f$ parameter is computed in a parameter dependent manner:
 * @f$ h = x \times e@f$, with @f$x \neq 0@f$ being the current parameter value.
 * If @f$x = 0@f$, @f$h = e@f$.
 * Default value is provided for @f$e@f$ and corresponds to the _h field.
 * The default value works fine in most cases, but you may want to change it using the setInterval method.
 *
 * @see AbstractNumericalDerivative, ThreePointsNumericalDerivative, FivePointsNumericalDerivative
 */
class TwoPointsNumericalDerivative :
  public AbstractNumericalDerivative
{
private:
  double f1_, f2_;

public:
  TwoPointsNumericalDerivative(Function* function) :
    AbstractNumericalDerivative(function),
    f1_(),
    f2_() {}
  TwoPointsNumericalDerivative(DerivableFirstOrder* function) :
    AbstractNumericalDerivative(function),
    f1_(),
    f2_() {}
  virtual ~TwoPointsNumericalDerivative() {}

  TwoPointsNumericalDerivative* clone() const { return new TwoPointsNumericalDerivative(*this); }

public:
  double getValue() const throw (Exception) { return f1_; }

  /**
   * @name The DerivableSecondOrder interface
   *
   * @{
   */
  double getSecondOrderDerivative(const std::string& variable) const
  throw (Exception)
  {
    throw Exception("Second order derivative not avalaible with two points method.");
  }

  double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const
  throw (Exception)
  {
    throw Exception("Unimplemented cross derivative.");
  }
  /** @} */

protected:
  void updateDerivatives(const ParameterList parameters)
  throw (ParameterNotFoundException, ConstraintException);
};
} // end of namespace bpp.

#endif // _TWOPOINTSNUMERICALDERIVATIVE_H_

