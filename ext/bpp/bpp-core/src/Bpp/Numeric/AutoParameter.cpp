//
// File: AutoParameter.cpp
// Created by: Julien Dutheil
// Created on: Tue Nov 11 22:15:16 2003
//

/*
 *  Copyright or © or Copr. CNRS, (November 17, 2004)
 * 
 *  This software is a computer program whose purpose is to provide classes
 *  for numerical calculus.
 * 
 *  This software is governed by the CeCILL  license under French law and
 *  abiding by the rules of distribution of free software.  You can  use, 
 *  modify and/ or redistribute the software under the terms of the CeCILL
 *  license as circulated by CEA, CNRS and INRIA at the following URL
 *  "http://www.cecill.info". 
 * 
 *  As a counterpart to the access to the source code and  rights to copy,
 *  modify and redistribute granted by the license, users are provided only
 *  with a limited warranty  and the software's author,  the holder of the
 *  economic rights,  and the successive licensors  have only  limited
 *  liability. 
 * 
 *  In this respect, the user's attention is drawn to the risks associated
 *  with loading,  using,  modifying and/or developing or reproducing the
 *  software by the user in light of its specific status of free software,
 *  that may mean  that it is complicated to manipulate,  and  that  also
 *  therefore means  that it is reserved for developers  and  experienced
 *  professionals having in-depth computer knowledge. Users are therefore
 *  encouraged to load and test the software's suitability as regards their
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and,  more generally, to use and operate it in the 
 *  same conditions as regards security. 
 * 
 *  The fact that you are presently reading this means that you have had
 *  knowledge of the CeCILL license and that you accept its terms.
 */

#include "AutoParameter.h"
#include "NumConstants.h"

#include <iostream>

using namespace std;

// Utils:
#include "../Text/TextTools.h"
#include "../App/ApplicationTools.h"

using namespace bpp;

/******************************************************************************/

string AutoParameter::CONSTRAINTS_AUTO   = "auto";
string AutoParameter::CONSTRAINTS_IGNORE = "ignore";
string AutoParameter::CONSTRAINTS_KEEP   = "keep";

/** Constructors: *************************************************************/

AutoParameter::AutoParameter(const std::string& name, double value, Constraint* constraint, bool attachConstraint)
noexcept(false):
Parameter(name, value, constraint, attachConstraint), messageHandler_(ApplicationTools::message) {}

AutoParameter::AutoParameter(const Parameter& p): Parameter(p), messageHandler_(ApplicationTools::message) {}

AutoParameter::AutoParameter(const AutoParameter& p): Parameter(p), messageHandler_(p.messageHandler_) {}

AutoParameter& AutoParameter::operator=(const AutoParameter& p)
{
  Parameter::operator=(p);
  messageHandler_ = p.messageHandler_;
  return *this;
}

/******************************************************************************/

void AutoParameter::setValue(double value) noexcept(false)
{
  try
  { 
    // First we try to assign this value:
    Parameter::setValue(value);
  }
  catch (ConstraintException& ce)
  { 
    // Aie, there's a pb here...
    if (messageHandler_)
    {
      (*messageHandler_) << "Constraint match at parameter ";
      (*messageHandler_) << name_;
      (*messageHandler_) << ", badValue = ";
      (*messageHandler_) << ce.getBadValue();
      (*messageHandler_) << " ";
      (*messageHandler_ << constraint_->getDescription()).endLine();
    }
    double limit = constraint_->getAcceptedLimit(value);
    try
    { // We try to assign the limit then.
      Parameter::setValue(limit);
    }
    catch(ConstraintException& ce2)
    { 
      // Aie, the limit is not reachable, so we perform a smaller step...
      //Parameter::setValue((getValue() + limit) / 2);
      try
      {
	// Try on the right:
	Parameter::setValue(limit + NumConstants::TINY());
      }
      catch(ConstraintException& ce3)
      {
	// Try on the left:
	Parameter::setValue(limit - NumConstants::TINY());
      }
    }
  }
}

/******************************************************************************/

