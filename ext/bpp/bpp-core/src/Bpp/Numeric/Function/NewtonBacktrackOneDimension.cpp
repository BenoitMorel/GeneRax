//
// File: NewtonBacktrackOneDimension.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 16 décembre 2010, à 15h 45
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

#include "NewtonBacktrackOneDimension.h"
#include "../NumTools.h"
#include "../../Text/TextTools.h"

using namespace bpp;

/******************************************************************************/

NewtonBacktrackOneDimension::NewtonBacktrackOneDimension(Function* function, double slope, double test) :
  AbstractOptimizer(function),
  fold_(0), f_(0), a_(0), alam_(0), alamin_(0), alam2_(0), b_(0), disc_(0), f2_(0), rhs1_(0), rhs2_(0), slope_(slope), test_(test), tmplam_(0)
    
{
  setDefaultStopCondition_(new NBODStopCondition(this));
  setStopCondition(*getDefaultStopCondition());
  setMaximumNumberOfEvaluations(10000);
}

/******************************************************************************/
  
void NewtonBacktrackOneDimension::doInit(const ParameterList& params) throw (Exception)
{
  // Set the initial value (no use here! Use setInitialValues() instead).
  if(params.size() != 1) throw Exception("NewtonBacktrackOneDimension::init(). This optimizer only deals with one parameter.");
  fold_ = getFunction()->f(getParameters());
  getStopCondition()->setTolerance(getStopCondition()->getTolerance()/test_);
  alamin_=getStopCondition()->getTolerance();
  alam_=1;
}

/******************************************************************************/

double NewtonBacktrackOneDimension::doStep() throw (Exception)
{
  if (alam_<alamin_){
    getParameter_(0).setValue(0);
    tolIsReached_=true;
    return fold_;
  }

  getParameter_(0).setValue(alam_);
  f_ = getFunction()->f(getParameters());

  if (f_<=fold_+alam_*0.0001*slope_){
    tolIsReached_=true;
    return f_;
  }

  if (alam_==1){
    tmplam_=-slope_/(2.0*(f_-fold_-slope_));
    f2_=f_;
    alam_=tmplam_>0.1?tmplam_:0.1;
    return f_;
  }
  
  rhs1_= f_-fold_-alam_*slope_;
  rhs2_= f2_-fold_-alam2_*slope_;

  a_=(rhs1_/(alam_*alam_)-rhs2_/(alam2_*alam2_))/(alam_-alam2_);
  b_=(-alam2_*rhs1_/(alam_*alam_)+alam_*rhs2_/(alam2_*alam2_))/(alam_-alam2_);

  if (a_==0.0)
    tmplam_= -slope_/(2.0*b_);
  else {
    disc_=b_*b_-3.0*a_*slope_;
    if (disc_<0.0)
      tmplam_=0.5*alam_;
    else
      if (b_<=0)
        tmplam_=(-b_+sqrt(disc_))/(3.0*a_);
      else
        tmplam_=-slope_/(b_+sqrt(disc_));
  }
  if (tmplam_> 0.5* alam_)
    tmplam_=0.5*alam_;

  alam2_=alam_;
  f2_=f_;
  alam_=tmplam_>0.1*alam_?tmplam_:0.1*alam_;

  return f_;
}

/******************************************************************************/

