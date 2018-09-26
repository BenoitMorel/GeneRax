//
// File: ParameterExceptions.h
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

#ifndef _PARAMETEREXCEPTIONS_H_
#define _PARAMETEREXCEPTIONS_H_

// From Utils:
#include "../Exceptions.h"

// From the STL:
#include <string>

namespace bpp
{

class Parameter;

/**
 * @brief The parameter exception base class.
 *
 * @see Exception
 */
class ParameterException :
  public Exception
{

	private:
		const Parameter* parameter_;
			
	public:	// Class constructors and destructor:
		/**
		 * @brief Build a new ParameterException object.
		 *
		 * @param text A message to be passed to the exception hierarchy.
		 * @param param A const pointer toward the parameter that threw the exception.
		 */	
		ParameterException(const std::string& text, const Parameter* param);

    ParameterException(const ParameterException& pe):
      Exception(pe),
      parameter_(pe.parameter_) {} //We explicitely do not want to make a hard copy of the pointer here.

    ParameterException& operator=(const ParameterException& pe)
    {
      Exception::operator=(pe);
      parameter_ = pe.parameter_;
      return *this;
    }
	
		virtual ~ParameterException() throw () {}
		
	public:
		/**
		 * @brief Get the parameter that threw the exception.
		 *
		 * @return A const pointer toward the parameter.
		 */
		virtual const Parameter * getParameter() const;
};

/**
 * @brief Exception thrown when a value do not match a given constraint.
 *
 * @see Constraint
 */
class ConstraintException : public ParameterException
{
	
	private:
		double badValue_;
	
	public: // Class constructor and destructor:
		/**
		 * @brief Build a new ConstraintException object.
		 *
		 * @param text     A message to be passed to the exception hierarchy.
		 * @param param    A const pointer toward the parameter that threw the exception.
		 * @param badValue The value that doesn't match the constraint.
		 */	
		ConstraintException(const std::string& text, const Parameter* param, double badValue);

		virtual ~ConstraintException() throw () {}
	
	public:
		/**
		 * @brief Get the value that doesn't match the constraint.
		 *
		 * @return The faulty value.
		 */
		virtual double getBadValue() const;
};

/*******************************************************************************/

/**
 * @brief Exception thrown when a parameter is not found,
 * for instance in a ParameterList.
 */
class ParameterNotFoundException : public Exception
{

	private:
		const std::string parameter_;
			
	public:	// Class constructors and destructor:
		/**
		 * @brief Build a new ParameterNotFoundException object.
		 *
		 * @param text     A message to be passed to the exception hierarchy.
		 * @param param    The name of the parameter not found.
		 */	
		ParameterNotFoundException(const std::string& text, const std::string& param = "");
	
		virtual ~ParameterNotFoundException() throw () {}
		
	public:
		/**
		 * @brief Get the name of the parameter not found.
		 *
		 * @return The parameter name.
		 */
		virtual std::string getParameter() const;		
};

} //end of namespace bpp.

#endif	//_PARAMETEREXCEPTIONS_H_

