//
// File: OptimizationStopCondition.cpp
// Created by: Julien Dutheil
// Created on: Tue Dec 23 11:51:31 2003
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 19, 2004)

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

#include "OptimizationStopCondition.h"
#include "Optimizer.h"
#include "../VectorTools.h"
#include "../NumTools.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

ParametersStopCondition::ParametersStopCondition(
  const Optimizer* optimizer) :
  AbstractOptimizationStopCondition(optimizer),
  lastParametersEstimates_(),
  newParametersEstimates_()
{
  init();
  if (newParametersEstimates_.size() == 0)
  {
    cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
         << "Be sure to have initialized the Optimizer first!" << endl; 
  }
}

ParametersStopCondition::ParametersStopCondition(
  const Optimizer* optimizer,
  double tolerance) :
  AbstractOptimizationStopCondition(optimizer, tolerance),
  lastParametersEstimates_(),
  newParametersEstimates_()
{
  init();
  if (newParametersEstimates_.size() == 0)
  {
    cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
         << "Be sure to have initialized the Optimizer first!" << endl; 
  }
}

ParametersStopCondition::ParametersStopCondition(
  const Optimizer* optimizer,
  int burnin) :
  AbstractOptimizationStopCondition(optimizer, burnin),
  lastParametersEstimates_(),
  newParametersEstimates_()
{
  init();
  if (newParametersEstimates_.size() == 0)
  {
    cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
         << "Be sure to have initialized the Optimizer first!" << endl; 
  }
}

ParametersStopCondition::ParametersStopCondition(
  const Optimizer* optimizer,
  double tolerance,
  int burnin) :
  AbstractOptimizationStopCondition(optimizer, tolerance, burnin),
  lastParametersEstimates_(),
  newParametersEstimates_()
{
  init();
  if (newParametersEstimates_.size() == 0)
  {
    cout << "DEBUG: WARNING!!! No parameter passed to ParametersStopCondition constructor. "
         << "Be sure to have initialized the Optimizer first!" << endl; 
  }
}

/******************************************************************************/

void ParametersStopCondition::init()
{
  AbstractOptimizationStopCondition::init();
  if (optimizer_->getFunction() != 0)
    newParametersEstimates_ = optimizer_->getParameters();
}

/******************************************************************************/

bool ParametersStopCondition::isToleranceReached() const
{
  callCount_++;
  lastParametersEstimates_ = newParametersEstimates_;
  newParametersEstimates_  = optimizer_->getParameters();
  if (callCount_ <= burnin_) return false;
  for (unsigned int i = 0; i < newParametersEstimates_.size(); i++)
  {
    Parameter& p = newParametersEstimates_[i];
    double lastEstimate = lastParametersEstimates_.getParameter(p.getName()).getValue();
    double newEstimate = p.getValue();
    double tol = NumTools::abs<double>(newEstimate - lastEstimate);
    if (tol > tolerance_)
    {
      return false;
    }
  }
  return true;
}

/******************************************************************************/

double ParametersStopCondition::getCurrentTolerance() const
{
  if (callCount_ > burnin_) {
    double maxTol = 0.;
    for (unsigned int i = 0; i < newParametersEstimates_.size(); i++)
    {
      Parameter& p = newParametersEstimates_[i];
      double lastEstimate = lastParametersEstimates_.getParameter(p.getName()).getValue();
      double newEstimate = p.getValue();
      double tol = NumTools::abs<double>(newEstimate - lastEstimate);
      if (tol > maxTol)
        maxTol = tol;
    }
    return maxTol;
  } else {
    return std::max(tolerance_, 1.);
  }
}

/******************************************************************************/

FunctionStopCondition::FunctionStopCondition(
  const Optimizer* optimizer) :
  AbstractOptimizationStopCondition(optimizer),
  lastFunctionValue_(-log(0.)),
  newFunctionValue_(-log(0.))
{
  init();
}

FunctionStopCondition::FunctionStopCondition(
  const Optimizer* optimizer,
  double tolerance) :
  AbstractOptimizationStopCondition(optimizer, tolerance),
  lastFunctionValue_(-log(0.)),
  newFunctionValue_(-log(0.))
{
  init();
}

FunctionStopCondition::FunctionStopCondition(
  const Optimizer* optimizer,
  int burnin) :
  AbstractOptimizationStopCondition(optimizer, burnin),
  lastFunctionValue_(-log(0.)),
  newFunctionValue_(-log(0.))
{
  init();
}

FunctionStopCondition::FunctionStopCondition(
  const Optimizer* optimizer,
  double tolerance,
  int burnin) :
  AbstractOptimizationStopCondition(optimizer, tolerance, burnin),
  lastFunctionValue_(-log(0.)),
  newFunctionValue_(-log(0.))
{
  init();
}

FunctionStopCondition::~FunctionStopCondition() {}

/******************************************************************************/

void FunctionStopCondition::init()
{
  AbstractOptimizationStopCondition::init();
  newFunctionValue_ = -log(0.);
  if (optimizer_->getFunction() != 0)
  {
    newFunctionValue_ = optimizer_->getFunctionValue();
  }
}

/******************************************************************************/

bool FunctionStopCondition::isToleranceReached() const
{
  callCount_++;
  lastFunctionValue_ = newFunctionValue_;
  newFunctionValue_  = optimizer_->getFunctionValue();
  if (callCount_ <= burnin_) return false;
  double tol = NumTools::abs<double>(newFunctionValue_ - lastFunctionValue_);
  return tol < tolerance_;
}

/******************************************************************************/

double FunctionStopCondition::getCurrentTolerance() const
{
  if (callCount_ > burnin_)
    return NumTools::abs<double>(newFunctionValue_ - lastFunctionValue_);
  else
    return std::max(tolerance_, 1.);
}

/******************************************************************************/

