//
// File: SimpleMultiDimensions.cpp
// Created by: Julien Dutheil
// Created on: Tue Nov 16 17:51 2004
//

/*
Copyright or © or Copr. CNRS, (November 19, 2004)

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

/******************************************************************************/

#include "SimpleMultiDimensions.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

SimpleMultiDimensions::SimpleMultiDimensions(Function* function):
  AbstractOptimizer(function), nbParams_(0), optimizer_(function)
{
  setDefaultStopCondition_(new FunctionStopCondition(this));
  setStopCondition(*getDefaultStopCondition());
  setOptimizationProgressCharacter("");
}

/******************************************************************************/

void SimpleMultiDimensions::setFunction(Function* function)
{
  AbstractOptimizer::setFunction(function);
  optimizer_.setFunction(function);
  getStopCondition()->init();
}

/******************************************************************************/

void SimpleMultiDimensions::doInit(const ParameterList& params) throw (Exception)
{
  nbParams_ = params.size();
  if (nbParams_ == 0) return;

  // Initialize optimizers:
  unsigned int nbEvalMax = nbEvalMax_ / static_cast<unsigned int>(nbParams_);
  optimizer_.setMaximumNumberOfEvaluations(nbEvalMax);
  optimizer_.setProfiler(getProfiler());
  optimizer_.setMessageHandler(getMessageHandler());
  optimizer_.getStopCondition()->setTolerance(getStopCondition()->getTolerance());
  optimizer_.setConstraintPolicy(getConstraintPolicy());
  optimizer_.setVerbose(getVerbose() > 0 ? getVerbose() - 1 : 0);
  optimizer_.setInitialInterval(0., 1.);
  getFunction()->setParameters(getParameters());
}

/******************************************************************************/

double SimpleMultiDimensions::doStep() throw (Exception)
{
  double f = getFunction()->getValue();
  for (unsigned int i = 0; i < nbParams_; i++)
  {
    if (getVerbose() > 0)
    {
      cout << getParameters()[i].getName() << ":";
      cout.flush();
    }
    // Re-init optimizer according to new values:
    double v = getParameters()[i].getValue();
    double t = std::max(0.000001, std::min(std::abs(v), getStopCondition()->getTolerance()));
    optimizer_.setInitialInterval(v - t, v + t);
    optimizer_.init(getParameters().subList(i));

    // Optimize through this dimension:
    f = optimizer_.optimize();
    if (getVerbose() > 0) cout << endl;
    getParameters_().matchParametersValues(getFunction()->getParameters());
    nbEval_ += optimizer_.getNumberOfEvaluations(); 
  }
  tolIsReached_ = nbParams_ <= 1;
  return f;
}

/******************************************************************************/

