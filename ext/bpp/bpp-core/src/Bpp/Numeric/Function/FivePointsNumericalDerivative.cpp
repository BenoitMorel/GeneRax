//
// File: FivePointsNumericalDerivative.cpp
// Created by: Julien Dutheil
// Created on: Thu Aug 17 15:00 2006
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

#include "FivePointsNumericalDerivative.h"
using namespace bpp;
using namespace std;

void FivePointsNumericalDerivative::updateDerivatives(const ParameterList parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  if (computeD1_ && variables_.size() > 0)
  {
    if (function1_)
      function1_->enableFirstOrderDerivatives(false);
    if (function2_)
      function2_->enableSecondOrderDerivatives(false);
    function_->setParameters(parameters);
    f3_ = function_->getValue();
    string lastVar;
    bool functionChanged = false;
    ParameterList p;
    bool start = true;
    for (unsigned int i = 0; i < variables_.size(); i++)
    {
      string var = variables_[i];
      if (!parameters.hasParameter(var))
        continue;
      if (!start)
      {
        vector<string> vars(2);
        vars[0] = var;
        vars[1] = lastVar;
        p = parameters.subList(vars);
      }
      else
      {
        p = parameters.subList(var);
        start = true;
      }
      lastVar = var;
      functionChanged = true;
      double value = function_->getParameterValue(var);
      double h = (1. + std::abs(value)) * h_;
      // Compute four other points:
      try
      {
        p[0].setValue(value - 2 * h);
        function_->setParameters(p);
        f1_ = function_->getValue();
        try
        {
          p[0].setValue(value + 2 * h);
          function_->setParameters(p);
          f5_ = function_->getValue();
          // No limit raised, use central approximation:
          p[0].setValue(value - h);
          function_->setParameters(p);
          f2_ = function_->getValue();
          p[0].setValue(value + h);
          function_->setParameters(p);
          f4_ = function_->getValue();
          der1_[i] = (f1_ - 8. * f2_ + 8. * f4_ - f5_) / (12. * h);
          der2_[i] = (-f1_ + 16. * f2_ - 30. * f3_ + 16. * f4_ - f5_) / (12. * h * h);
        }
        catch (ConstraintException& ce)
        {
          // Right limit raised, use backward approximation:
          p[0].setValue(value - h);
          function_->setParameters(p);
          f2_ = function_->getValue();
          p[0].setValue(value - 2 * h);
          function_->setParameters(p);
          f1_ = function_->getValue();
          der1_[i] = (f3_ - f2_) / h;
          der2_[i] = (f3_ - 2. * f2_ + f1_) / (h * h);
        }
      }
      catch (ConstraintException& ce)
      {
        // Left limit raised, use forward approximation:
        p[0].setValue(value + h);
        function_->setParameters(p);
        f4_ = function_->getValue();
        p[0].setValue(value + 2 * h);
        function_->setParameters(p);
        f5_ = function_->getValue();
        der1_[i] = (f4_ - f3_) / h;
        der2_[i] = (f5_ - 2. * f4_ + f3_) / (h * h);
      }
    }
    // Reset last parameter and compute analytical derivatives if any.
    if (function1_)
      function1_->enableFirstOrderDerivatives(computeD1_);
    if (function2_)
      function2_->enableSecondOrderDerivatives(computeD2_);
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
    function_->setParameters(parameters);
    // Just in  case derivatives are not computed:
    f3_ = function_->getValue();
  }
}

