//
// File: BrentOneDimension.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov 17 11:45:58 2003
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#include "BrentOneDimension.h"
#include "OneDimensionOptimizationTools.h"
#include "../NumTools.h"
#include "../NumConstants.h"
#include "../../Text/TextTools.h"

using namespace bpp;

/******************************************************************************/

bool BrentOneDimension::BODStopCondition::isToleranceReached() const
{
  callCount_++;
  if (callCount_ <= burnin_) return false;
  const BrentOneDimension* bod = dynamic_cast<const BrentOneDimension *>(optimizer_);
  return getCurrentTolerance() <= (bod->tol2 - 0.5 * (bod->b - bod->a));
}
    
/******************************************************************************/

double BrentOneDimension::BODStopCondition::getCurrentTolerance() const
{
  // NRC Test for done:
  const BrentOneDimension* bod = dynamic_cast<const BrentOneDimension *>(optimizer_);
  return NumTools::abs(bod->x - bod->xm);
}
 
/******************************************************************************/

BrentOneDimension::BrentOneDimension(Function* function) :
  AbstractOptimizer(function),
  a(0), b(0), d(0), e(0), etemp(0), fu(0), fv(0), fw(0), fx(0), p(0), q(0), r(0), tol1(0), tol2(0),
  u(0), v(0), w(0), x(0), xm(0), _xinf(0), _xsup(0), isInitialIntervalSet_(false)
{
  setDefaultStopCondition_(new BODStopCondition(this));
  setStopCondition(*getDefaultStopCondition());
  setMaximumNumberOfEvaluations(10000);
}

/******************************************************************************/
  
double BrentOneDimension::ZEPS  = 1.0e-10;
  
/******************************************************************************/
  
void BrentOneDimension::doInit(const ParameterList& params) throw (Exception)
{
  if (params.size() != 1)
    throw Exception("BrentOneDimension::init(). This optimizer only deals with one parameter.");

  // Bracket the minimum.
  Bracket bracket = OneDimensionOptimizationTools::bracketMinimum(_xinf, _xsup, getFunction(), getParameters());
  if (getVerbose() > 0)
  {
    printMessage("Initial bracketing:");
    printMessage("A: x = " + TextTools::toString(bracket.a.x) + ", f = " + TextTools::toString(bracket.a.f));
    printMessage("B: x = " + TextTools::toString(bracket.b.x) + ", f = " + TextTools::toString(bracket.b.f));
    printMessage("C: x = " + TextTools::toString(bracket.c.x) + ", f = " + TextTools::toString(bracket.c.f));
  }
  
  // This will be the distance moved on the step before last.
  e = 0.0;

  // a and b must be in ascending order, but input abscissa need not be.
  a = (bracket.a.x < bracket.c.x ? bracket.a.x : bracket.c.x);
  b = (bracket.a.x > bracket.c.x ? bracket.a.x : bracket.c.x);
  // Initializations...
  fw = fv = fx = getFunction()->f(getParameters());
  if (fx < bracket.b.f)
  {
    //We don't want to lose our initial guess!
    x = w = v = bracket.b.x = getParameters()[0].getValue();
  }
  else
  {
    x = w = v = bracket.b.x;
    getParameter_(0).setValue(x);
    fw = fv = fx = getFunction()->f(getParameters());
  }
}

/******************************************************************************/
  
void BrentOneDimension::setInitialInterval(double inf, double sup)
{
  if(sup > inf)
  {
    _xinf = inf; _xsup = sup;
  }
  else
  {
    _xinf = sup; _xsup = inf;
  }
  isInitialIntervalSet_ = true;
}

/******************************************************************************/

double BrentOneDimension::doStep() throw (Exception)
{
  xm   = 0.5 * (a + b);
  tol2 = 2.0 * (tol1 = getStopCondition()->getTolerance() * NumTools::abs(x) + ZEPS);
  
  if(NumTools::abs(e) > tol1)
  {
    r = (x - w) * (fx - fv);
    q = (x - v) * (fx - fw);
    p = (x - v) * q - (x - w) * r;
    q = 2.0 * (q - r);
    if (q > 0.0) p = -p;
    q = NumTools::abs(q);
    etemp = e;
    e = d;
    if (NumTools::abs(p) >= NumTools::abs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x))
      d = NumConstants::GOLDEN_RATIO_C() * (e = (x >= xm ? a - x : b - x));
    else
    {
      d = p / q;
      u = x + d;
      if (u - a < tol2 || b - u < tol2)
        d = NumTools::sign(tol1, xm - x);
    }
  }
  else
  {
    d = NumConstants::GOLDEN_RATIO_C() * (e = (x >= xm ? a - x : b - x));
  }
  u = (NumTools::abs(d) >= tol1 ? x + d : x + NumTools::sign(tol1, d));

  // Function evaluaton:
  ParameterList pl = getParameters();
  pl[0].setValue(u);

  fu = getFunction()->f(pl);

  if (fu <= fx)
  {
    if (u >= x) a = x; else b = x;
    NumTools::shift(v, w, x, u);
    NumTools::shift(fv, fw, fx, fu);
  }
  else
  {
    if (u < x) a = u; else b = u;
    if (fu <= fw || w == x)
    {
      v  = w;
      w  = u;
      fv = fw;
      fw = fu;
    }
    else if (fu <= fv || v == x || v == w)
    {
      v  = u;
      fv = fu;
    }
  }

  // Store results for this step:
  getParameter_(0).setValue(x);
  return fx;
}

/******************************************************************************/
  
double BrentOneDimension::optimize() throw (Exception)
{
  if (!isInitialIntervalSet_)
    throw Exception("BrentOneDimension::optimize. Initial interval not set: call the 'setInitialInterval' method first!");
  AbstractOptimizer::optimize();
  // Apply parameters and evaluate function at minimum point:
  currentValue_ = getFunction()->f(getParameters());
  return currentValue_;
}

/******************************************************************************/

