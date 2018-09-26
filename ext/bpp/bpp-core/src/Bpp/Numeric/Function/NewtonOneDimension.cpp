//
// File: NewtonOneDimension.cpp
// Created by: Julien Dutheil
// Created on: Thu Apr 26 14:16 2007
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

#include "NewtonOneDimension.h"
#include "../NumTools.h"
#include "../../Text/TextTools.h"

using namespace bpp;

/******************************************************************************/

NewtonOneDimension::NewtonOneDimension(DerivableSecondOrder* function) :
  AbstractOptimizer(function),
  _param(),
  _maxCorrection(10)
{
  setDefaultStopCondition_(new FunctionStopCondition(this));
  setStopCondition(*getDefaultStopCondition());
  nbEvalMax_ = 10000;
}

/******************************************************************************/
  
void NewtonOneDimension::doInit(const ParameterList& params) throw (Exception)
{
  // Set the initial value (no use here! Use setInitialValues() instead).
  if (params.size() != 1)
    throw Exception("NewtonOneDimension::init(). This optimizer only deals with one parameter.");
  _param = params[0].getName();
  currentValue_ = getFunction()->f(getParameters());  
  getStopCondition()->init();
}

/******************************************************************************/

double NewtonOneDimension::doStep() throw (Exception)
{
  double movement;  
  ParameterList newPoint = getParameters();
  ParameterList bckPoint = getFunction()->getParameters();
  double newValue;
  double  firstOrderDerivative = getFunction()->getFirstOrderDerivative(_param);
  double secondOrderDerivative = getFunction()->getSecondOrderDerivative(_param);
  if (secondOrderDerivative <= 0)
  {
    printMessage("!!! Second order derivative is negative (" + TextTools::toString(getParameters()[0].getValue()) + "). No move performed.");
    //movements[i] = 0;  // We want to reach a minimum, not a maximum!
    // My personnal improvement:
    movement = -firstOrderDerivative / secondOrderDerivative;
  }
  else movement = firstOrderDerivative / secondOrderDerivative;
  if (std::isnan(movement))
  {
    printMessage("!!! Non derivable point. No move performed. (f=" + TextTools::toString(currentValue_) + ", d1=" + TextTools::toString(firstOrderDerivative) + ", d2=" + TextTools::toString(secondOrderDerivative) + ").");
    movement = 0; // Either first or second order derivative is infinity. This may happen when the function == inf at this point.
  }
  newPoint[0].setValue(getParameters()[0].getValue() - movement); 
  newValue = getFunction()->f(newPoint);

  // Check newValue:
  unsigned int count = 0;
  while (newValue > currentValue_)
  {
    //Restore previous point (all parameters in case of global constraint):
    getFunction()->setParameters(bckPoint);

    count++;
    if (count >= _maxCorrection)
    {
      printMessage("!!! Felsenstein-Churchill correction applied too much time. Stopping here. Convergence probably not reached.");
      tolIsReached_ = true;
      return currentValue_;
      //throw Exception("NewtonOneDimension::step(). Felsenstein-Churchill correction applied more than 10000 times.");
    }
    printMessage("!!! Function at new point is greater than at current point: " + TextTools::toString(newValue) + ">" + TextTools::toString(currentValue_) + ". Applying Felsenstein-Churchill correction, value = " + TextTools::toString(newPoint[0].getValue()));
    movement = movement / 2;
    newPoint[0].setValue(getParameters()[0].getValue() - movement);
    newValue = getFunction()->f(newPoint);
  }

  getParameters_() = newPoint; // Function as been set to newPoint by the call of f(newPoint).
  return newValue;
}

/******************************************************************************/

