//
// File: AbstractOptimizer.cpp
// Created by: Julien Dutheil
// Created on: Mon Dec 22 12:18:09 2003
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

#include "AbstractOptimizer.h"
#include "../AutoParameter.h"
#include "../../Text/TextTools.h"
#include "../../App/ApplicationTools.h"

// From the STL:
#include <iomanip>
#include <time.h>

using namespace std;
using namespace bpp;

/******************************************************************************/

AbstractOptimizer::AbstractOptimizer(Function* function):
  function_(function),
  parameters_(),
  messageHandler_(ApplicationTools::message),
  profiler_(ApplicationTools::message),
  constraintPolicy_(AutoParameter::CONSTRAINTS_KEEP),
  stopCondition_(0), defaultStopCondition_(0),
  verbose_(true), isInitialized_(false), startTime_(), listeners_(),
  updateParameters_(false), stepChar_("*"),
  nbEvalMax_(1000000), nbEval_(0),
  currentValue_(0), tolIsReached_(false)
{
}

/******************************************************************************/

AbstractOptimizer::AbstractOptimizer(const AbstractOptimizer& opt):
  function_(opt.function_),
  parameters_(opt.parameters_),
  messageHandler_(opt.messageHandler_),
  profiler_(opt.profiler_),
  constraintPolicy_(opt.constraintPolicy_),
  stopCondition_(0), defaultStopCondition_(0),
  verbose_(opt.verbose_),
  isInitialized_(opt.isInitialized_),
  startTime_(opt.startTime_),
  listeners_(), //We do not copy listeners!
  updateParameters_(opt.updateParameters_),
  stepChar_(opt.stepChar_),
  nbEvalMax_(opt.nbEvalMax_),
  nbEval_(opt.nbEval_),
  currentValue_(opt.currentValue_),
  tolIsReached_(opt.tolIsReached_)
{
  if (opt.stopCondition_)
    {
      stopCondition_        = dynamic_cast<OptimizationStopCondition *>(opt.stopCondition_->clone());
      stopCondition_->setOptimizer(this);
    }
  else
    stopCondition_        = 0;
  if (opt.defaultStopCondition_)
    {
      defaultStopCondition_ = dynamic_cast<OptimizationStopCondition *>(opt.defaultStopCondition_->clone());
      defaultStopCondition_->setOptimizer(this);
    }
  else
    defaultStopCondition_ = 0;
  //In case of AutoParameter instances, we must actualize the pointers toward _messageHandler:
  if (isInitialized_)
    {
      if(constraintPolicy_ == AutoParameter::CONSTRAINTS_AUTO)   autoParameter();
      else if(constraintPolicy_ == AutoParameter::CONSTRAINTS_IGNORE) ignoreConstraints();
    }
}

/******************************************************************************/

AbstractOptimizer& AbstractOptimizer::operator=(const AbstractOptimizer& opt)
{
  function_               = opt.function_;
  parameters_             = opt.parameters_;
  messageHandler_         = opt.messageHandler_;
  profiler_               = opt.profiler_;
  constraintPolicy_       = opt.constraintPolicy_;
  tolIsReached_           = opt.tolIsReached_;
  if (opt.stopCondition_)
  {
    stopCondition_        = dynamic_cast<OptimizationStopCondition *>(opt.stopCondition_->clone());
    stopCondition_->setOptimizer(this);
  }
  else
    stopCondition_        = 0;
  if (opt.defaultStopCondition_)
  {
    defaultStopCondition_ = dynamic_cast<OptimizationStopCondition *>(opt.defaultStopCondition_->clone());
    defaultStopCondition_->setOptimizer(this);
  }
  else
    defaultStopCondition_ = 0;
  nbEvalMax_              = opt.nbEvalMax_;
  nbEval_                 = opt.nbEval_;
  verbose_                = opt.verbose_;
  isInitialized_          = opt.isInitialized_;
  //In case of AutoParameter instances, we must actualize the pointers toward messageHandler_:
  if (isInitialized_)
  {
    if (constraintPolicy_ == AutoParameter::CONSTRAINTS_AUTO)   autoParameter();
    else if (constraintPolicy_ == AutoParameter::CONSTRAINTS_IGNORE) ignoreConstraints();
  }
  startTime_              = opt.startTime_;
  listeners_.resize(0); //Reset listener list, do not copy it!
  updateParameters_       = opt.updateParameters_;
  stepChar_               = opt.stepChar_;
  return *this;
}

/******************************************************************************/
	
void AbstractOptimizer::init(const ParameterList& params) throw (Exception)
{
  if (!function_) throw Exception("AbstractOptimizer::init. Optimizer currently has no function.");
  //We do this in order to keep original constraints:
  parameters_ = params;
  //More secure, but too slow:
  //parameters_ = function_->getParameters().subList(params.getParameterNames());
  //parameters_.matchParametersValues(params);
  if (constraintPolicy_ == AutoParameter::CONSTRAINTS_AUTO) autoParameter();
  else if (constraintPolicy_ == AutoParameter::CONSTRAINTS_IGNORE) ignoreConstraints();
  doInit(params);
  nbEval_ = 0;
  tolIsReached_ = false;
  isInitialized_ = true;
  time(&startTime_);
  currentValue_ = function_->getValue();

  profile("Step\t");
  for (unsigned int i = 0; i < parameters_.size(); i++)
  {
    profile(parameters_[i].getName() + "\t"); 
  }
  profileln("Function\tTime");

  //Parameters must be assigned by doInit:

  printPoint(parameters_, currentValue_);
  
  // Initialize the StopCondition:
  stopCondition_->init();
  fireOptimizationInitializationPerformed(OptimizationEvent(this));
}

/******************************************************************************/

double AbstractOptimizer::step() throw (Exception)
{
  currentValue_ = doStep();
  printPoint(parameters_, currentValue_);
  fireOptimizationStepPerformed(OptimizationEvent(this));
  if (listenerModifiesParameters())
  {
    if (!updateParameters_)
      parameters_.matchParametersValues(function_->getParameters());
    //else already done!
     //_currentValue = function_->getValue();
    //Often useless, but avoid some bizare behaviour in particular cases:
    currentValue_ = function_->f(parameters_);
  }
  tolIsReached_ = tolIsReached_ || stopCondition_->isToleranceReached();
  return currentValue_;
}

/**************************************************************************/

double AbstractOptimizer::optimize() throw (Exception)
{
  if (!isInitialized_)
    throw Exception("AbstractOptimizer::optimize. Optimizer not initialized: call the 'init' method first!");
  tolIsReached_ = false;
  for (nbEval_ = 1; nbEval_ < nbEvalMax_ && ! tolIsReached_; nbEval_++)
  {
    if (verbose_)
      ApplicationTools::displayUnlimitedGauge(nbEval_, "Optimizing... ");
    step();
  }
  return currentValue_;
}

/******************************************************************************/

void AbstractOptimizer::profile(double v)
{
  if (profiler_) *profiler_ << v;
}	

/******************************************************************************/

void AbstractOptimizer::profileln(double v)
{
  if (profiler_) (*profiler_ << v).endLine();
}
	
/******************************************************************************/

void AbstractOptimizer::profile(unsigned int v)
{
  if (profiler_) *profiler_ << v;
}
/******************************************************************************/

void AbstractOptimizer::profileln(unsigned int v)
{
  if (profiler_) (*profiler_ << v).endLine();
}
	
/******************************************************************************/

void AbstractOptimizer::profile(const std::string& s)
{
  if (profiler_) *profiler_ << s;
}	

/******************************************************************************/

void AbstractOptimizer::profileln(const std::string& s)
{
  if (profiler_) (*profiler_ << s).endLine();
}
	
/******************************************************************************/

void AbstractOptimizer::printPoint(const ParameterList& params, double value)
{
  size_t ndim = params.size();
  profile(nbEval_);
  profile("\t");
  for (size_t j = 0; j < ndim; j++)
  {
    profile(TextTools::toString(params[j].getValue()));
    profile("\t"); 
  }
  profile(value);
  profile("\t");
  time_t seconds;
  time(&seconds);
  profileln(difftime(seconds, startTime_));
}

/******************************************************************************/

void AbstractOptimizer::printMessage(const std::string& message)
{
  if (messageHandler_) (*messageHandler_ << message).endLine();
}

/******************************************************************************/

void AbstractOptimizer::autoParameter()
{
  for (unsigned int i = 0; i < parameters_.size(); i++)
  {
    AutoParameter ap(parameters_[i]);
    ap.setMessageHandler(messageHandler_);
    parameters_.setParameter(i, ap);
  }
}

/******************************************************************************/

void AbstractOptimizer::ignoreConstraints()
{
  for (unsigned int i = 0; i < parameters_.size(); i++)
  {
    parameters_[i].removeConstraint();
  }
}

/******************************************************************************/

void AbstractOptimizer::fireOptimizationInitializationPerformed(const OptimizationEvent& event)
{
  for (unsigned int i = 0; i < listeners_.size(); i++)
  {
    listeners_[i]->optimizationInitializationPerformed(event);
  }
}

/******************************************************************************/

void AbstractOptimizer::fireOptimizationStepPerformed(const OptimizationEvent& event)
{
  for (unsigned int i = 0; i < listeners_.size(); i++)
  {
    listeners_[i]->optimizationStepPerformed(event);
  }
}

/******************************************************************************/

bool AbstractOptimizer::listenerModifiesParameters() const
{
  for (unsigned int i = 0; i < listeners_.size(); i++)
  {
    if (listeners_[i]->listenerModifiesParameters())
      return true;
  }
  return false;
}

/******************************************************************************/

