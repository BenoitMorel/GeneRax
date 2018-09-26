//
// File: MetaOptimizer.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 12 16:05 2007
// From file: NewtonBrentMetaOptimizer.cpp
// Created on: Tue Nov 17 17:22 2004
// 
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

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

/**************************************************************************/

#include "MetaOptimizer.h"
#include "../../App/ApplicationTools.h"

using namespace bpp;
using namespace std;

/**************************************************************************/

string MetaOptimizerInfos::IT_TYPE_STEP = "step";
string MetaOptimizerInfos::IT_TYPE_FULL = "full";

/**************************************************************************/

MetaOptimizer::MetaOptimizer(
    Function* function,
    MetaOptimizerInfos* desc,
    unsigned int n):
  AbstractOptimizer(function),
  optDesc_(desc), optParameters_(desc->getNumberOfOptimizers()),
  nbParameters_(desc->getNumberOfOptimizers()), n_(n),
  precisionStep_(-1.), stepCount_(0), initialValue_(-1.)
{
  setDefaultStopCondition_(new FunctionStopCondition(this));
  setStopCondition(*getDefaultStopCondition());
  precisionStep_ = log10(getStopCondition()->getTolerance()) / n_;
  setOptimizationProgressCharacter("");
}

/**************************************************************************/

MetaOptimizer::MetaOptimizer(
    const MetaOptimizer& opt):
  AbstractOptimizer(opt),
  optDesc_(dynamic_cast<MetaOptimizerInfos *>(opt.optDesc_->clone())),
  optParameters_(opt.optParameters_),
  nbParameters_(opt.nbParameters_),
  n_(opt.n_),
  precisionStep_(opt.precisionStep_),
  stepCount_(opt.stepCount_),
  initialValue_(opt.initialValue_)
{}

/**************************************************************************/

MetaOptimizer& MetaOptimizer::operator=(
    const MetaOptimizer& opt)
{
  AbstractOptimizer::operator=(opt);
  optDesc_       = dynamic_cast<MetaOptimizerInfos *>(opt.optDesc_->clone());
  optParameters_ = opt.optParameters_;
  nbParameters_  = opt.nbParameters_;
  n_             = opt.n_;
  precisionStep_ = opt.precisionStep_;
  stepCount_     = opt.stepCount_;
  initialValue_  = opt.initialValue_;
  return *this;
}

/**************************************************************************/

MetaOptimizer::~MetaOptimizer()
{
  // Delete all optimizers:
  delete optDesc_;
}

/**************************************************************************/

void MetaOptimizer::doInit(const ParameterList& parameters)
  throw (Exception)
{
  optParameters_.resize(optDesc_->getNumberOfOptimizers());
  for (unsigned int i = 0; i < optDesc_->getNumberOfOptimizers(); i++)
  {
    optParameters_[i].reset();
    for (size_t j = 0; j < optDesc_->getParameterNames(i).size(); j++)
    {
      string pname = optDesc_->getParameterNames(i)[j];
      if (parameters.hasParameter(pname))
      {
        optParameters_[i].addParameter(parameters.getParameter(pname));
      }
    }
    nbParameters_[i] = optParameters_[i].size();
  }

  // Initialize optimizers:
  for (unsigned int i = 0; i < optDesc_->getNumberOfOptimizers(); i++)
  {
    if (nbParameters_[i] > 0)
    {
      Optimizer * opt = optDesc_->getOptimizer(i);
      dynamic_cast<AbstractOptimizer*>(opt)->updateParameters(updateParameters());
      opt->setProfiler(getProfiler());
      opt->setMessageHandler(getMessageHandler());
      opt->setConstraintPolicy(getConstraintPolicy());
      opt->setVerbose(getVerbose() > 0 ? getVerbose() - 1 : 0);
    }
  }
  
  // Actualize parameters:
  getParameters_().matchParametersValues(getFunction()->getParameters());
  
  getFunction()->setParameters(getParameters());
  initialValue_ = getFunction()->getValue();
  // Reset counter:
  stepCount_ = 1;
  // Recompute step if precision has changed:
  precisionStep_ = (log10(getStopCondition()->getTolerance()) - log10(initialValue_)) / n_;
}

/**************************************************************************/

double MetaOptimizer::doStep() throw (Exception)
{
  stepCount_++;
  
  int tolTest = 0;
  double tol = getStopCondition()->getTolerance();
  if (stepCount_ <= n_)
  {
    tol = initialValue_ * pow(10, stepCount_ * precisionStep_);
  }
  
  for (unsigned int i = 0; i < optDesc_->getNumberOfOptimizers(); i++)
  {
    if (nbParameters_[i] > 0)
    {
      if (getVerbose() > 1 && ApplicationTools::message)
      {
        (ApplicationTools::message->endLine() << optDesc_->getName(i)).endLine();
        ApplicationTools::message->flush();
      }
      if (optDesc_->requiresFirstOrderDerivatives(i))
        dynamic_cast<DerivableFirstOrder*>(getFunction())->enableFirstOrderDerivatives(true);
      if (optDesc_->requiresSecondOrderDerivatives(i))  
        dynamic_cast<DerivableSecondOrder*>(getFunction())->enableSecondOrderDerivatives(true);

      optParameters_[i].matchParametersValues(getParameters());
      Optimizer * opt = optDesc_->getOptimizer(i);
      opt->getStopCondition()->setTolerance(tol);
      opt->init(optParameters_[i]);
      if (optDesc_->getIterationType(i) == MetaOptimizerInfos::IT_TYPE_STEP)
        opt->step();
      else if (optDesc_->getIterationType(i) == MetaOptimizerInfos::IT_TYPE_FULL)
        opt->optimize();
      else throw Exception("MetaOptimizer::step. Unknown iteration type specified.");
      nbEval_ += opt->getNumberOfEvaluations();
      if (optDesc_->requiresFirstOrderDerivatives(i))
        dynamic_cast<DerivableFirstOrder*>(getFunction())->enableFirstOrderDerivatives(false);
      if (optDesc_->requiresSecondOrderDerivatives(i))  
        dynamic_cast<DerivableSecondOrder*>(getFunction())->enableSecondOrderDerivatives(false);
      if (getVerbose() > 1) cout << endl;
      getParameters_().matchParametersValues(opt->getParameters());
    }
    tolTest += nbParameters_[i] > 0 ? 1 : 0;
  }
  tolIsReached_ = (tolTest == 1);
   
  return getFunction()->getValue();
}

/**************************************************************************/

