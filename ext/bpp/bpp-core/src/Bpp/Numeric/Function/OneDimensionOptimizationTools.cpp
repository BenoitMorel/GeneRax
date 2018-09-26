//
// File: OneDimensionOptimizationTools.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov 17 11:15:22 2003
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for numerical calculus. This file is part of the Bio++ project.

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

#include "NewtonBacktrackOneDimension.h"
#include "BrentOneDimension.h"
#include "OneDimensionOptimizationTools.h"
#include "../NumTools.h"
#include "../NumConstants.h"

using namespace bpp;
using namespace std;

/******************************************************************************
*                              The Point class                               *
******************************************************************************/
inline void BracketPoint::set(double xval, double fval) { this->x = xval; this->f = fval; }

/******************************************************************************
*                             The Bracket class                              *
******************************************************************************/
inline void Bracket::setA(double xa, double fa) { a.set(xa, fa); }
inline void Bracket::setB(double xb, double fb) { b.set(xb, fb); }
inline void Bracket::setC(double xc, double fc) { c.set(xc, fc); }

/******************************************************************************/

Bracket OneDimensionOptimizationTools::bracketMinimum(
  double a,
  double b,
  Function* function,
  ParameterList parameters)
{
  Bracket bracket;
  // Copy the parameter to use.
  bracket.a.x = a;
  parameters[0].setValue(bracket.a.x); bracket.a.f = function->f(parameters);
  bracket.b.x = b;
  parameters[0].setValue(bracket.b.x); bracket.b.f = function->f(parameters);

  while (std::isnan(bracket.b.f)|| std::isinf(bracket.b.f))
  {
    bracket.b.x /= 1.1;
    parameters[0].setValue(bracket.b.x); bracket.b.f = function->f(parameters);
  }
  
  if (bracket.b.f > bracket.a.f)
  {
    // Switch roles of first and second point so that we can go downhill
    // in the direction from a to b.
    NumTools::swap<double>(bracket.a.x, bracket.b.x);
    NumTools::swap<double>(bracket.a.f, bracket.b.f);
  }

  // First guess for third point:
  bracket.c.x = bracket.b.x + NumConstants::GOLDEN_RATIO_PHI() * (bracket.b.x - bracket.a.x);
  parameters[0].setValue(bracket.c.x); bracket.c.f = function->f(parameters);


  // Keep returning here until we bracket:
  while (bracket.b.f > bracket.c.f)
  {
    // Compute xu by parabolic extrapolation from a, b, c. TINY is used to prevent
    // any possible division by 0.
    double r = (bracket.b.x - bracket.a.x) * (bracket.b.f - bracket.c.f);
    double q = (bracket.b.x - bracket.c.x) * (bracket.b.f - bracket.a.f);

    double xu = bracket.b.x - ((bracket.b.x - bracket.c.x) * q - (bracket.b.x - bracket.a.x) * r) /
                (2.0 * NumTools::sign(NumTools::max(NumTools::abs(q - r), NumConstants::VERY_TINY()), q - r));
    double xulim = (bracket.b.x) + GLIMIT * (bracket.c.x - bracket.b.x);
    double fu;

    // We don't go farther than this.
    // Test various possibilities:
    if ((bracket.b.x - xu) * (xu - bracket.c.x) > 0.0)
    {
      parameters[0].setValue(xu); fu = function->f(parameters);
      if (fu < bracket.c.f)
      {
        bracket.setA(bracket.b.x, bracket.b.f);
        bracket.setB(xu, fu);
        return bracket;
      }
      else if (fu > bracket.b.f)
      {
        bracket.setC(xu, fu);
        return bracket;
      }
      // Parabolic fit was no use.
      // Use default magnification.
      xu = bracket.c.x + NumConstants::GOLDEN_RATIO_PHI() * (bracket.c.x - bracket.b.x);
      parameters[0].setValue(xu); fu = function->f(parameters);
    }
    else if ((bracket.c.x - xu) * (xu - xulim) > 0.0)
    {
      // Parabolic fit is between point 3 and its allowed limit.
      parameters[0].setValue(xu); fu = function->f(parameters);
      if (fu < bracket.c.f)
      {
        NumTools::shift<double>(bracket.b.x, bracket.c.x, xu, bracket.c.x + NumConstants::GOLDEN_RATIO_PHI() * (bracket.c.x - bracket.b.x));
        parameters[0].setValue(xu);
        NumTools::shift<double>(bracket.b.f, bracket.c.f, fu, function->f(parameters));
      }
    }
    else if ((xu - xulim) * (xulim - bracket.c.x) >= 0.0)
    {
      // Limit parabolic xu to maximum allowed value.
      xu = xulim;
      parameters[0].setValue(xu); fu = function->f(parameters);
    }
    else
    {
      // Reject parabolic xu, use default magnification.
      xu = bracket.c.x + NumConstants::GOLDEN_RATIO_PHI() * (bracket.c.x - bracket.b.x);
      parameters[0].setValue(xu); fu = function->f(parameters);
    }
    // Eliminate oldest point and continue.
    NumTools::shift<double>(bracket.a.x, bracket.b.x, bracket.c.x, xu);
    NumTools::shift<double>(bracket.a.f, bracket.b.f, bracket.c.f, fu);
  }
  return bracket;
}

/******************************************************************************/

unsigned int OneDimensionOptimizationTools::lineMinimization(
  DirectionFunction& f1dim,
  ParameterList& parameters,
  std::vector<double>& xi,
  double tolerance,
  OutputStream* profiler,
  OutputStream* messenger,
  unsigned int verbose)
{
  // Initial guess for brackets:
  double ax = 0.;
  double xx = 0.01;

  f1dim.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  f1dim.setMessageHandler(messenger);
  f1dim.init(parameters, xi);
  BrentOneDimension bod(&f1dim);
  bod.setMessageHandler(messenger);
  bod.setProfiler(profiler);
  bod.setVerbose(verbose >= 1 ? 1 : 0);
  bod.setOptimizationProgressCharacter(".");
  bod.getStopCondition()->setTolerance(0.01);
  bod.setInitialInterval(ax, xx);
  bod.setConstraintPolicy(AutoParameter::CONSTRAINTS_KEEP);
  ParameterList singleParameter;
  singleParameter.addParameter(Parameter("x", 0.0));
  bod.init(singleParameter);
  bod.optimize();
  // Update parameters:
  // parameters.matchParametersValues(f1dim.getFunction()->getParameters());

  double xmin = f1dim.getParameters()[0].getValue();
  for (size_t j = 0; j < parameters.size(); j++)
  {
    xi[j] *= xmin;
    parameters[j].setValue(parameters[j].getValue() + xi[j]);
  }
  return bod.getNumberOfEvaluations();
}

/******************************************************************************/

unsigned int OneDimensionOptimizationTools::lineSearch(DirectionFunction& f1dim,
                                                       ParameterList& parameters,
                                                       std::vector<double>& xi,
                                                       std::vector<double>& gradient,
                                                       OutputStream* profiler,
                                                       OutputStream* messenger,
                                                       unsigned int verbose)
{
  size_t size = xi.size();

  f1dim.setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
  f1dim.setMessageHandler(messenger);
  f1dim.init(parameters, xi);

  double slope = 0;
  for (size_t i = 0; i < size; i++)
  {
    slope += xi[i] * gradient[i];
  }

  //  if (slope>=0)
  //  throw Exception("Slope problem in OneDimensionOptimizationTools::lineSearch. Slope="+TextTools::toString(slope));

  double x, temp, test = 0;
  for (size_t i = 0; i < size; ++i)
  {
    x = abs(parameters[i].getValue());
    temp = abs(xi[i]);
    if (x > 1.0)
      temp /= x;
    if (temp > test)
      test = temp;
  }

  NewtonBacktrackOneDimension nbod(&f1dim, slope, test);

  nbod.setMessageHandler(messenger);
  nbod.setProfiler(profiler);
  nbod.setVerbose(verbose >= 1 ? 1 : 0);
  nbod.setOptimizationProgressCharacter(".");
  nbod.getStopCondition()->setTolerance(0.0001);
  //  nbod.setInitialInterval(ax, xx);
  nbod.setConstraintPolicy(AutoParameter::CONSTRAINTS_KEEP);
  ParameterList singleParameter;
  singleParameter.addParameter(Parameter("x", 0.0));
  nbod.init(singleParameter);
  nbod.optimize();
  // Update parameters:
  // parameters.matchParametersValues(f1dim.getFunction()->getParameters());

  double xmin = f1dim.getParameters()[0].getValue();
  for (unsigned int j = 0; j < parameters.size(); j++)
  {
    xi[j] *= xmin;
    parameters[j].setValue(parameters[j].getValue() + xi[j]);
  }

  return nbod.getNumberOfEvaluations();
}

/******************************************************************************/

double OneDimensionOptimizationTools::GLIMIT = 100.0;

/******************************************************************************/
