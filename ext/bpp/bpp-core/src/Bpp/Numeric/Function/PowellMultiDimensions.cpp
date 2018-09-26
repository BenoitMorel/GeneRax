//
// File: PowellMultiDimensions.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov 17 15:16:45 2003
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

#include "PowellMultiDimensions.h"
#include "BrentOneDimension.h"
#include "OneDimensionOptimizationTools.h"
#include "../NumTools.h"

using namespace bpp;

/******************************************************************************/

bool PowellMultiDimensions::PMDStopCondition::isToleranceReached() const {
  callCount_++;
  if (callCount_ <= burnin_) return false;
  return getCurrentTolerance() < tolerance_;
}

double PowellMultiDimensions::PMDStopCondition::getCurrentTolerance() const
{
  // NRC Test for done:
  const PowellMultiDimensions* pmd = dynamic_cast<const PowellMultiDimensions*>(optimizer_);
  double fp   = pmd->fp_;
  double fret = pmd->fret_;
  return 2.0 * NumTools::abs(fp - fret) / (NumTools::abs(fp) + NumTools::abs(fret));
}
    
/******************************************************************************/

PowellMultiDimensions::PowellMultiDimensions(Function* function) :
AbstractOptimizer(function), fp_(0), fret_(0), pt_(), xi_(), ncom_(0), pcom_(), xicom_(), f1dim_(function)
{
  setDefaultStopCondition_(new PMDStopCondition(this));
  setStopCondition(*getDefaultStopCondition());
}

/******************************************************************************/

void PowellMultiDimensions::doInit(const ParameterList& params) throw (Exception)
{
  // Build the initial matrix:
  size_t n = params.size();
  xi_.resize(n);
  for (size_t i = 0; i < n; i++)
  {
    // Copy the parameter list:
    xi_[i].resize(n);
    for(unsigned int j = 0; j < n; j++)
    {
      // Set the directions to unit vectors:
      xi_[i][j] = (j == i) ? 1 : 0;
    }
  }
  
  // Starting point:
  fret_ = getFunction()->f(getParameters());
  pt_   = getParameters();
}
  
/******************************************************************************/
  
double PowellMultiDimensions::doStep() throw (Exception)
{
  size_t n = getParameters().size();
  fp_ = fret_;
  unsigned int ibig = 0;
  double del = 0.0; // Will be the biggest function decrease
  Vdouble xit(n);
  
  // In each iteration, loop over all directions in the set.
  double fptt;
  for(unsigned int i = 0; i < n; i++)
  {
    // Copy the direction:
    for(unsigned int j = 0; j < n; j++)
    {
      xit[j] = xi_[j][i];
    }
    fptt = fret_;
    nbEval_ += OneDimensionOptimizationTools::lineMinimization(f1dim_,
        getParameters_(), xit, getStopCondition()->getTolerance(),
        0, getMessageHandler(), getVerbose() > 0 ? getVerbose() - 1 : 0);
    fret_ = getFunction()->f(getParameters());
    if (getVerbose() > 2) printPoint(getParameters(), fret_);
    if (fret_ > fp_) throw Exception("DEBUG: PowellMultiDimensions::doStep(). Line minimization failed!");
    if (fptt - fret_ > del)
    {
      del = fptt - fret_;
      ibig = i;
    }
  }

  ParameterList ptt = getParameters();
  for (unsigned int j = 0; j < n; j++)
  {
    ptt[j].setValue(2.0 * getParameters()[j].getValue() - pt_[j].getValue());
    xit[j] = getParameters()[j].getValue() - pt_[j].getValue();
    pt_[j].setValue(getParameters()[j].getValue());
  }
  fptt = getFunction()->f(ptt);
  if (fptt < fp_)
  {
    double t = 2.0 * (fp_ - 2.0 * fret_ + fptt) * NumTools::sqr(fp_ - fret_ - del) - del * NumTools::sqr(fp_ - fptt);
    if (t < 0.0)
    {
      //cout << endl << "New direction: drection " << ibig << " removed." << endl;
      nbEval_ += OneDimensionOptimizationTools::lineMinimization(f1dim_,
          getParameters_(), xit, getStopCondition()->getTolerance(),
          0, getMessageHandler(), getVerbose() > 0 ? getVerbose() - 1 : 0);
      fret_ = getFunction()->f(getParameters());
      if (fret_ > fp_) throw Exception("DEBUG: PowellMultiDimensions::doStep(). Line minimization failed!");
      for (unsigned int j = 0; j < n; j++)
      {
        xi_[j][ibig]  = xi_[j][n - 1];
        xi_[j][n - 1] = xit[j];
      }
    }
  }
  else getFunction()->setParameters(getParameters());

  return fret_;
}

/******************************************************************************/

double PowellMultiDimensions::optimize() throw (Exception)
{
  AbstractOptimizer::optimize();
  // Apply best parameter:
  return getFunction()->f(getParameters());
}

/******************************************************************************/

