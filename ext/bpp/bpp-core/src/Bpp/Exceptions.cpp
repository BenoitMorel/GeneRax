//
// File Exceptions.cpp
// Created by: Guillaume Deuchst
//              Julien Dutheil
//              Sylvain Gaillard
// Last modification : Thu Jul 22 2004
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide utilitary
classes. This file belongs to the Bio++ Project.

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

#include "Exceptions.h"
#include "Text/TextTools.h"

using namespace bpp;
using namespace std;

/******************************************************************************/
  
BadIntegerException::BadIntegerException(const char* text, int badInt):
  Exception(string(text) + "(" + TextTools::toString(badInt) + ")"),
  badInt_(badInt) {}

BadIntegerException::BadIntegerException(const std::string& text, int badInt):
  Exception(text + "(" + TextTools::toString(badInt) + ")"),
  badInt_(badInt) {}

/******************************************************************************/

BadNumberException::BadNumberException(const char* text, double badNumber):
  Exception(string(text) + "(" + TextTools::toString(badNumber) + ")"),
  badNumber_(badNumber) {}

BadNumberException::BadNumberException(const std::string& text, double badNumber):
  Exception(text + "(" + TextTools::toString(badNumber) + ")"),
  badNumber_(badNumber) {}

/******************************************************************************/

NumberFormatException::NumberFormatException(const char* text, const std::string& badNumber):
  Exception(string(text) + "(" + badNumber + ")"),
  badNumber_(badNumber) {}

NumberFormatException::NumberFormatException(const std::string& text, const std::string& badNumber):
  Exception(text + "(" + badNumber + ")"),
  badNumber_(badNumber) {}

/******************************************************************************/

vector<size_t> IndexOutOfBoundsException::getBounds() const
{
  vector<size_t> bounds(2);
  bounds[0] = lowerBound_;
  bounds[1] = upperBound_;
  return bounds;
}

/******************************************************************************/
  
BadSizeException::BadSizeException(const std::string& text, size_t badSize, size_t correctSize):
  Exception("Incorrect size " + TextTools::toString(badSize) + ", expected " + TextTools::toString(correctSize) + ". " + text),
  badSize_(badSize),
  correctSize_(correctSize) {}

/******************************************************************************/

OutOfRangeException::OutOfRangeException(const std::string& text, double badValue, double lowerBound, double upperBound):
  Exception(TextTools::toString(badValue) + " out of [" + TextTools::toString(lowerBound) + ", " + TextTools::toString(upperBound) +  "])" + text),
  lowerBound_(lowerBound),
  upperBound_(upperBound)
{}

/******************************************************************************/

