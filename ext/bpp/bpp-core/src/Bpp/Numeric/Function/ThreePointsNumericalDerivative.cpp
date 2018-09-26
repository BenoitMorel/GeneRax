//
// File: ThreePointsNumericalDerivative.cpp
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

#include "ThreePointsNumericalDerivative.h"

using namespace bpp;
using namespace std;

void ThreePointsNumericalDerivative::updateDerivatives(const ParameterList parameters)
throw (ParameterNotFoundException, ConstraintException)
{
  if (computeD1_ && variables_.size() > 0)
  {
    if (function1_)
      function1_->enableFirstOrderDerivatives(false);
    if (function2_)
      function2_->enableSecondOrderDerivatives(false);
    function_->setParameters(parameters);
    f2_ = function_->getValue();
    if ((abs(f2_) >= NumConstants::VERY_BIG()) || std::isnan(f2_))
    {
      for (size_t i = 0; i < variables_.size(); ++i)
      {
        der1_[i] = log(-1);
        der2_[i] = log(-1);
      }
      return;
    }

    string lastVar;
    bool functionChanged = false;
    ParameterList p;
    bool start = true;
    for (size_t i = 0; i < variables_.size(); ++i)
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
        start = false;
      }
      lastVar = var;
      functionChanged = true;
      double value = function_->getParameterValue(var);
      double h = -(1. + std::abs(value)) * h_;
      if (abs(h) < p[0].getPrecision())
        h = h < 0 ? -p[0].getPrecision() : p[0].getPrecision();
      double hf1(0), hf3(0);
      unsigned int nbtry = 0;

      // Compute f1_
      while (hf1 == 0)
      {
        try
        {
          p[0].setValue(value + h);
          function_->setParameters(p); // also reset previous parameter...

          p = p.subList(0);
          f1_ = function_->getValue();
          if ((abs(f1_) >= NumConstants::VERY_BIG()) || std::isnan(f1_))
            throw ConstraintException("f1_ too large", &p[0], f1_);
          else
            hf1 = h;
        }
        catch (ConstraintException& ce)
        {
          if (++nbtry == 10) // no possibility to compute derivatives
            break;
          else if (h < 0)
            h = -h;  // try on the right
          else
            h /= -2;  // try again on the left with smaller interval
        }
      }

      if (hf1 != 0)
      {
        // Compute f3_
        if (h < 0)
          h = -h;  // on the right
        else
          h /= 2;  //  on the left with smaller interval

        nbtry = 0;
        while (hf3 == 0)
        {
          try
          {
            p[0].setValue(value + h);
            function_->setParameters(p); // also reset previous parameter...

            p = p.subList(0);
            f3_ = function_->getValue();
            if ((abs(f3_) >= NumConstants::VERY_BIG()) || std::isnan(f3_))
              throw ConstraintException("f3_ too large", &p[0], f3_);
            else
              hf3 = h;
          }
          catch (ConstraintException& ce)
          {
            if (++nbtry == 10) // no possibility to compute derivatives
              break;
            else if (h < 0)
              h = -h;  // try on the right
            else
              h /= -2;  // try again on the left with smaller interval
          }
        }
      }

      if (hf3 == 0)
      {
        der1_[i] = log(-1);
        der2_[i] = log(-1);
      }
      else
      {
        der1_[i] = (f1_ - f3_) / (hf1 - hf3);
        der2_[i] = ((f1_ - f2_) / hf1 - (f3_ - f2_) / hf3) * 2 / (hf1 - hf3);
      }
    }


    if (computeCrossD2_)
    {
      string lastVar1, lastVar2;
      for (unsigned int i = 0; i < variables_.size(); i++)
      {
        string var1 = variables_[i];
        if (!parameters.hasParameter(var1))
          continue;
        for (unsigned int j = 0; j < variables_.size(); j++)
        {
          if (j == i)
          {
            crossDer2_(i, j) = der2_[i];
            continue;
          }
          string var2 = variables_[j];
          if (!parameters.hasParameter(var2))
            continue;

          vector<string> vars(2);
          vars[0] = var1;
          vars[1] = var2;
          if (i > 0 && j > 0)
          {
            if (lastVar1 != var1 && lastVar1 != var2)
              vars.push_back(lastVar1);
            if (lastVar2 != var1 && lastVar2 != var2)
              vars.push_back(lastVar2);
          }
          p = parameters.subList(vars);

          double value1 = function_->getParameterValue(var1);
          double value2 = function_->getParameterValue(var2);
          double h1 = (1. + std::abs(value1)) * h_;
          double h2 = (1. + std::abs(value2)) * h_;

          // Compute 4 additional points:
          try
          {
            p[0].setValue(value1 - h1);
            p[1].setValue(value2 - h2);
            function_->setParameters(p); // also reset previous parameter...
            vector<size_t> tmp(2);
            tmp[0] = 0;
            tmp[1] = 1;
            p = p.subList(tmp); // removed the previous parameters.
            f11_ = function_->getValue();

            p[1].setValue(value2 + h2);
            function_->setParameters(p.subList(1));
            f12_ = function_->getValue();

            p[0].setValue(value1 + h1);
            function_->setParameters(p.subList(0));
            f22_ = function_->getValue();

            p[1].setValue(value2 - h2);
            function_->setParameters(p.subList(1));
            f21_ = function_->getValue();

            crossDer2_(i, j) = ((f22_ - f21_) - (f12_ - f11_)) / (4 * h1 * h2);
          }
          catch (ConstraintException& ce)
          {
            throw Exception("ThreePointsNumericalDerivative::setParameters. Could not compute cross derivatives at limit.");
          }

          lastVar1 = var1;
          lastVar2 = var2;
        }
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
    // Just in case derivatives are not computed:
    f2_ = function_->getValue();
  }
}

