//
// File: ParameterExceptions.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov  3 18:05:36 2003
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

#include "ParameterExceptions.h"
#include "Parameter.h"

// From Utils:
#include "../Text/TextTools.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

ParameterException::ParameterException(const std::string& text, const Parameter* param) :
	Exception("ParameterException: " + text + (param != 0 ? "(" + param -> getName() + ")" : string(""))),
	parameter_(param) {};
		
const Parameter* ParameterException::getParameter() const { return parameter_; }

/******************************************************************************/

ConstraintException::ConstraintException(const std::string& text, const Parameter* param, double badValue) :
	ParameterException("ConstraintException: " + text + "(" + TextTools::toString(badValue) + ")"
      + (param->hasConstraint() ? param->getConstraint()->getDescription() : "no constraint"), param),
	badValue_(badValue) {};
		
double ConstraintException::getBadValue() const { return badValue_; }

/******************************************************************************/

ParameterNotFoundException::ParameterNotFoundException(const string& text, const string& param) :
	Exception("ParameterNotFoundException: " + text + (!TextTools::isEmpty(param) ? "(" + param + ")" : string(""))),
	parameter_(param) {};
		
std::string ParameterNotFoundException::getParameter() const { return parameter_; }

/******************************************************************************/

