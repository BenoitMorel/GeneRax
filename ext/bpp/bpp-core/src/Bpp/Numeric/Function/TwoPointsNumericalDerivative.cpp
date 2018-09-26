//
// File: TwoPointsNumericalDerivative.cpp
// Created by: Julien Dutheil
// Created on: Mon May 28 10:33 2007
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#include "TwoPointsNumericalDerivative.h"

using namespace bpp;
using namespace std;

void TwoPointsNumericalDerivative::updateDerivatives(const ParameterList parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  if (computeD1_ && variables_.size() > 0)
  {
    if (function1_)
      function1_->enableFirstOrderDerivatives(false);
    if (function2_)
      function2_->enableSecondOrderDerivatives(false);
    function_->setParameters(parameters);
    f1_ = function_->getValue();
    string lastVar;
    bool functionChanged = false;
    bool start = true;
    for (unsigned int i = 0; i < variables_.size(); i++)
    {
      string var = variables_[i];
      if (!parameters.hasParameter(var))
        continue;
      ParameterList p;
      if (!start)
      {
        vector<string> vars(2);
        vars[0] = var;
        vars[1] = lastVar;
        lastVar = var;
        functionChanged = true;
        p = parameters.subList(vars);
      }
      else
      {
        p = parameters.subList(var);
        lastVar = var;
        functionChanged = true;
        start = false;
      }
      double value = function_->getParameterValue(var);
      double h = (1 + std::abs(value)) * h_;
      // Compute one other point:
      try
      {
        p[0].setValue(value + h);
        function_->setParameters(p);
        f2_ = function_->getValue();
      }
      catch (ConstraintException& ce1)
      {
        // Right limit raised, use backward approximation:
        try
        {
          p[0].setValue(value - h);
          function_->setParameters(p);
          f2_ = function_->getValue();
          der1_[i] = (f1_ - f2_) / h;
        }
        catch (ConstraintException& ce2)
        {
          // PB: can't compute derivative, because of a two narrow interval (lower than h)
          throw ce2;
        }
      }
      // No limit raised, use forward approximation:
      der1_[i] = (f2_ - f1_) / h;
    }
    // Reset last parameter and compute analytical derivatives if any:
    if (function1_)
      function1_->enableFirstOrderDerivatives(computeD1_);
    if (functionChanged)
      function_->setParameters(parameters.subList(lastVar));
  }
  else
  {
    // Reset initial value and compute analytical derivatives if any.
    if (function1_)
      function1_->enableFirstOrderDerivatives(computeD1_);
    if (function2_)
      function2_->enableSecondOrderDerivatives(computeD2_);
    // Just in case derivatives are not computed:
    function_->setParameters(parameters);
    f1_ = function_->getValue();
  }
}

