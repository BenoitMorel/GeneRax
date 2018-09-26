//
// File: Parameter.cpp
// Created by: Julien Dutheil
// Created on: Wed Oct 15 15:40:47 2003
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

#include "Parameter.h"
#include <cmath>

//From Utils:
#include "../Text/TextTools.h"

using namespace bpp;

#include <iostream>
using namespace std;

/******************************************************************************/

ParameterEvent::ParameterEvent(Parameter* parameter): parameter_(parameter) {}

/** Constructors: *************************************************************/

Parameter::Parameter(const std::string& name, double value, Constraint* constraint, bool attachConstraint, double precision)
noexcept(false) :
  name_(name), value_(0), precision_(0), constraint_(constraint), attach_(attachConstraint), listeners_(), listenerAttach_()
{
  // This may throw a ConstraintException:
  setValue(value);
  setPrecision(precision);
}

Parameter::Parameter(const std::string& name, double value, const Constraint* constraint, double precision)
noexcept(false) :
  name_(name), value_(0), precision_(0), constraint_(constraint ? constraint->clone() : 0), attach_(true), listeners_(), listenerAttach_()
{
  // This may throw a ConstraintException:
  setValue(value);
  setPrecision(precision);
}

Parameter::Parameter(const Parameter& p) :
  name_(p.name_),
  value_(p.value_),
  precision_(p.precision_),
  constraint_(0),
  attach_(p.attach_),
  listeners_(p.listeners_),
  listenerAttach_(p.listenerAttach_)
{
  if (p.attach_ && p.constraint_)
    constraint_ = p.constraint_->clone();
  else
    constraint_ = p.constraint_;
  for (size_t i = 0; i < listeners_.size(); i++)
    if (listenerAttach_[i])
      listeners_[i] = dynamic_cast<ParameterListener*>(p.listeners_[i]->clone());
}

Parameter& Parameter::operator=(const Parameter& p)
{
  name_           = p.name_;
  value_          = p.value_;
  precision_      = p.precision_;
  attach_         = p.attach_;
  if (p.attach_ && p.constraint_)
    constraint_   = p.constraint_->clone();
  else
    constraint_   = p.constraint_;
  listeners_      = p.listeners_;
  listenerAttach_ = p.listenerAttach_;
  for (size_t i = 0; i < listeners_.size(); i++)
    if (listenerAttach_[i])
      listeners_[i] = dynamic_cast<ParameterListener*>(p.listeners_[i]->clone());
  return *this;	
}

/** Destructor: ***************************************************************/

Parameter::~Parameter()
{
  if (attach_ && constraint_) delete constraint_;
  for (size_t i = 0; i < listeners_.size(); i++)
    if (listenerAttach_[i])
      delete listeners_[i];
} 

/** Value: ********************************************************************/

void Parameter::setValue(double value) noexcept(false)
{
  if (std::abs(value-value_)>precision_/2){
    if (constraint_ && !constraint_->isCorrect(value)) 
      throw ConstraintException("Parameter::setValue", this, value);
    value_ = value;
    ParameterEvent event(this);
    fireParameterValueChanged(event);
  }
}

/** Precision: ********************************************************************/

void Parameter::setPrecision(double precision)
{
  precision_=(precision<0)?0:precision;
}

/** Constraint: ***************************************************************/

void Parameter::setConstraint(Constraint* constraint, bool attach)
{
  if (!constraint){
    if (constraint_ && attach_)
      delete constraint_;
    constraint_=0;
    attach_=false;
  }
  else {
    if (constraint->isCorrect(value_)) {
      if (constraint_ && attach_)
        delete constraint_;
      constraint_ = constraint;
      attach_= attach;
    }
  else throw ConstraintException("Parameter::setConstraint", this, value_);
  }
}


const Constraint* Parameter::removeConstraint()
{
  const Constraint * c = constraint_;
  constraint_ = 0;
  return c;
}

/******************************************************************************/

void Parameter::removeParameterListener(const std::string& listenerId)
{
  for (unsigned int i = 0; i < listeners_.size(); i++)
    {
      if (listeners_[i]->getId() == listenerId)
        {
          if (listenerAttach_[i]) delete listeners_[i];
          listeners_.erase(listeners_.begin() + i);
          listenerAttach_.erase(listenerAttach_.begin() + i);
        }
    }
}

/******************************************************************************/

bool Parameter::hasParameterListener(const std::string& listenerId)
{
  for (unsigned int i = 0; i < listeners_.size(); i++)
    if (listeners_[i]->getId() == listenerId)
      return true;
  return false;
}

/******************************************************************************/

const IntervalConstraint Parameter::R_PLUS(true, 0, true);

const IntervalConstraint Parameter::R_PLUS_STAR(true, 0, false);

const IntervalConstraint Parameter::R_MINUS(false, 0, true);

const IntervalConstraint Parameter::R_MINUS_STAR(false, 0, false);

const IntervalConstraint Parameter::PROP_CONSTRAINT_IN(0, 1, true, true);

const IntervalConstraint Parameter::PROP_CONSTRAINT_EX(0, 1, false, false);

/******************************************************************************/

