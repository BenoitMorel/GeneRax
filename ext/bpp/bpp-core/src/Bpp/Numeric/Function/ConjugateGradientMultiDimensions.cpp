//
// File: ConjugateGradientMultiDimensions.cpp
// Created by: Julien Dutheil
// Created on: Wed Apr 11 16:51 2007
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

#include "ConjugateGradientMultiDimensions.h"
#include "OneDimensionOptimizationTools.h"

#include "../VectorTools.h"

using namespace bpp;

/******************************************************************************/

ConjugateGradientMultiDimensions::ConjugateGradientMultiDimensions(DerivableFirstOrder* function):
  AbstractOptimizer(function), optimizer_(function),
  xi_(), h_(), g_(), f1dim_(function)
{
  setDefaultStopCondition_(new FunctionStopCondition(this));
  setStopCondition(*getDefaultStopCondition());
}

/******************************************************************************/

void ConjugateGradientMultiDimensions::doInit(const ParameterList & params) throw (Exception)
{
  size_t nbParams = params.size();
  g_.resize(nbParams);
  h_.resize(nbParams);
  xi_.resize(nbParams);
  getFunction_()->enableFirstOrderDerivatives(true);
  getFunction_()->setParameters(params);
  getGradient(xi_);
  for(size_t i = 0; i < nbParams; i++)
  {
    g_[i]  = -xi_[i];
    xi_[i] = h_[i] = g_[i];
  }
}

/******************************************************************************/

double ConjugateGradientMultiDimensions::doStep() throw (Exception)
{
  double gg, gam, f, dgg;
  size_t n = getParameters().size();
  //Loop over iterations.
  getFunction_()->enableFirstOrderDerivatives(false);
  nbEval_ += OneDimensionOptimizationTools::lineMinimization(f1dim_,
      getParameters_(), xi_, getStopCondition()->getTolerance(),
      0, 0, getVerbose() > 0 ? getVerbose() - 1 : 0);

  getFunction_()->enableFirstOrderDerivatives(true);
  f = getFunction()->f(getParameters());

  if (tolIsReached_)
  {
    return f;
  }
  getGradient(xi_);

  dgg = gg = 0.0;
  for (unsigned j = 0; j < n; j++)
  {
    gg += g_[j] * g_[j];
    /* dgg += xi[j] * xi[j]; */ //This statement for Fletcher-Reeves.
    dgg += (xi_[j] + g_[j]) * xi_[j]; //This statement for Polak-Ribiere.
  }
  
  if (gg == 0.0)
  { 
    //Unlikely. If gradient is exactly zero then
    return f;
  }
  gam = dgg / gg;

  if (!(std::isnan(gam) || std::isinf(gam)))
  {
    for(unsigned int j = 0; j < n; j++)
    {
      g_[j] = -xi_[j];
      xi_[j] = h_[j] = g_[j] + gam * h_[j];
    }
  }
  
  return f;
}

/******************************************************************************/

void ConjugateGradientMultiDimensions::getGradient(std::vector<double>& gradient) const
{
  for (size_t i = 0; i < gradient.size(); ++i)
  {
    gradient[i] = getFunction()->getFirstOrderDerivative(getParameters()[i].getName());
  }
}

/******************************************************************************/

