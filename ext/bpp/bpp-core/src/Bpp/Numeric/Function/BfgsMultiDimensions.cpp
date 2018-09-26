//
// File: BfgsMultiDimensions.cpp
// Created by: Laurent Guéguen
// Created on: Dec 16 13:49 2010
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

#include "BfgsMultiDimensions.h"
#include "OneDimensionOptimizationTools.h"

using namespace bpp;

/******************************************************************************/

BfgsMultiDimensions::BfgsMultiDimensions(DerivableFirstOrder* function) :
  AbstractOptimizer(function),
  //gtol_(gtol),
  slope_(0),
  Up_(),
  Lo_(),
  p_(),
  gradient_(),
  xi_(),
  dg_(),
  hdg_(),
  hessian_(),
  f1dim_(function)
{
  setDefaultStopCondition_(new FunctionStopCondition(this));
  setStopCondition(*getDefaultStopCondition());
  setOptimizationProgressCharacter(".");
}

/******************************************************************************/

void BfgsMultiDimensions::doInit(const ParameterList& params) throw (Exception)
{
  size_t nbParams = params.size();
  p_.resize(nbParams);
  gradient_.resize(nbParams);
  xi_.resize(nbParams);
  dg_.resize(nbParams);
  hdg_.resize(nbParams);
  Up_.resize(nbParams);
  Lo_.resize(nbParams);

  hessian_.resize(nbParams);
  for (size_t i = 0; i < nbParams; i++)
  {
    hessian_[i].resize(nbParams);
  }

  for (size_t i = 0; i < nbParams; i++)
  {
    const Constraint* cp = params[i].getConstraint();
    if (!cp)
    {
      Up_[i] = NumConstants::VERY_BIG();
      Lo_[i] = -NumConstants::VERY_BIG();
    }
    else
    {
      Up_[i] = cp->getAcceptedLimit(NumConstants::VERY_BIG()) - NumConstants::TINY();
      Lo_[i] = cp->getAcceptedLimit(-NumConstants::VERY_BIG()) + NumConstants::TINY();
    }
  }

  getFunction_()->enableFirstOrderDerivatives(true);
  getFunction_()->setParameters(params);

  getGradient(gradient_);

  for (size_t i = 0; i < nbParams; i++)
  {
    p_[i] = getParameters()[i].getValue();

    for (unsigned int j = 0; j < nbParams; j++)
    {
      hessian_[i][j] = 0.0;
    }
    hessian_[i][i] = 1.0;
  }


  double sum = 0;
  for (size_t i = 0; i < nbParams; i++)
  {
    sum += p_[i] * p_[i];
  }
}

/******************************************************************************/

double BfgsMultiDimensions::doStep() throw (Exception)
{
  double f;
  size_t n = getParameters().size();
  // Loop over iterations.

  unsigned int i;

  for (i = 0; i < n; i++)
  {
    p_[i] = getParameters()[i].getValue();
  }

  setDirection();

  getFunction()->enableFirstOrderDerivatives(false);
  nbEval_ += OneDimensionOptimizationTools::lineSearch(f1dim_,
                                                       getParameters_(), xi_,
                                                       gradient_,
                                                       // getStopCondition()->getTolerance(),
                                                       0, 0,
                                                       getVerbose() > 0 ? getVerbose() - 1 : 0);
  getFunction()->enableFirstOrderDerivatives(true);

  for (i = 0; i < n; i++)
  {
    xi_[i] = getParameters_()[i].getValue() - p_[i];
  }

  f = getFunction()->f(getParameters());
  if (f > currentValue_) {
    printMessage("!!! Function increase !!!");
    printMessage("!!! Optimization might have failed. Try to reparametrize your function to remove constraints.");
    tolIsReached_ = true;
    return f;
  }

  if (tolIsReached_)
  {
    return f;
  }

  //double temp, test = 0.0;
  //for (i = 0; i < n; i++)
  //{
  //  temp = xi_[i];
  //  if (p_[i] > 1.0)
  //    temp /= p_[i];
  //  if (temp > test)
  //    test = temp;
  //}

  //if (test < 1e-7)
  //{
  //  tolIsReached_ = true;
  //  return f;
  //}

  for (i = 0; i < n; i++)
  {
    dg_[i] = gradient_[i];
  }

  getGradient(gradient_);
  //test = 0.0;

  //for (i = 0; i < n; i++)
  //{
  //  temp = abs(gradient_[i]);
  //  if (abs(p_[i]) > 1.0)
  //    temp /= p_[i];
  //  if (temp > test)
  //    test = temp;
  //}

  //if (f > 1.0)
  //  test /= f;

  //if (test < gtol_)
  //{
  //  tolIsReached_ = true;
  //  return f;
  //}

  for (i = 0; i < n; i++)
  {
    dg_[i] = gradient_[i] - dg_[i];
  }

  for (i = 0; i < n; i++)
  {
    hdg_[i] = 0;
    for (unsigned int j = 0; j < n; j++)
    {
      hdg_[i] += hessian_[i][j] * dg_[j];
    }
  }

  double fae(0), fac(0), sumdg(0), sumxi(0);

  for (i = 0; i < n; i++)
  {
    fac += dg_[i] * xi_[i];
    fae += dg_[i] * hdg_[i];
    sumdg += dg_[i] * dg_[i];
    sumxi += xi_[i] * xi_[i];
  }

  if (fac > sqrt(1e-7 * sumdg * sumxi))
  {
    fac = 1.0 / fac;
    double fad = 1.0 / fae;
    for (i = 0; i < n; i++)
    {
      dg_[i] = fac * xi_[i] - fad * hdg_[i];
    }
    for (i = 0; i < n; i++)
    {
      for (unsigned int j = i; j < n; j++)
      {
        hessian_[i][j] += fac * xi_[i] * xi_[j] - fad * hdg_[i] * hdg_[j] + fae * dg_[i] * dg_[j];
        hessian_[j][i] = hessian_[i][j];
      }
    }
  }

  return f;
}

/******************************************************************************/

void BfgsMultiDimensions::getGradient(std::vector<double>& gradient) const
{
  for (unsigned int i = 0; i < gradient.size(); i++)
  {
    gradient[i] = getFunction()->getFirstOrderDerivative(getParameters()[i].getName());
  }
}

/******************************************************************************/

void BfgsMultiDimensions::setDirection()
{
  size_t nbParams = getParameters().size();

  for (size_t i = 0; i < nbParams; i++)
  {
    xi_[i] = 0;
    for (unsigned int j = 0; j < nbParams; j++)
    {
      xi_[i] -= hessian_[i][j] * gradient_[j];
    }
  }

  double v = 1, alpmax = 1;
  for (size_t i = 0; i < nbParams; i++)
  {
    if ((xi_[i] > 0) && (p_[i] + NumConstants::TINY() * xi_[i] < Up_[i]))
      v = (Up_[i] - p_[i]) / xi_[i];
    else if ((xi_[i] < 0) && (p_[i] + NumConstants::TINY() * xi_[i] > Lo_[i]))
      v = (Lo_[i] - p_[i]) / xi_[i];
    if (v < alpmax)
      alpmax = v;
  }

  for (size_t i = 0; i < nbParams; i++)
  {
    if (p_[i] + NumConstants::TINY() * xi_[i] >= Up_[i])
      xi_[i] = Up_[i] - p_[i];
    else if (p_[i] + NumConstants::TINY() * xi_[i] <= Lo_[i])
      xi_[i] = Lo_[i] - p_[i];
    else
      xi_[i] *= alpmax;
  }
}
