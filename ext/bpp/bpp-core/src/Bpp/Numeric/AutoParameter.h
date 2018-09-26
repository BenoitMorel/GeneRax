//
// File: AutoParameter.h
// Created by: Julien Dutheil
// Created on: Tue Nov 11 22:15:16 2003
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

#ifndef _AUTOPARAMETER_H_
#define _AUTOPARAMETER_H_

#include "Parameter.h"

//From Utils:
#include "../Io/OutputStream.h"

namespace bpp
{

/**
 * @brief The AutoParameter class.
 *
 * This class overides the setValue() method of class Parameter so that no
 * Exception is thrown. This allows to perform optimization under constraint.
 */ 

class AutoParameter:
  public Parameter
{
	private:
    OutputStream* messageHandler_;
	
	public:
		
		/**
		 * @brief Build a new AutoParameter.
		 *
		 * @param name The parameter name.
		 * @param value The parameter value.
		 * @param constraint An optional pointer toward a Constraint object.
     * @param attachConstraint Tell if the constraint must be attached to this parameter, or shared
     * between different objects. See Parameter.
		 * @throw ConstraintException If the parameter value does not match the contraint.
		 */
		AutoParameter(const std::string& name = "", double value = 0, Constraint* constraint = 0, bool attachConstraint = false)
      noexcept(false);

		/**
		 * @brief Copy constructor.
		 *
		 * @param param The parameter to copy.
		 */
		AutoParameter(const Parameter& param);
	
		/**
		 * @brief Copy constructor.
		 *
		 * @param param The parameter to copy.
		 */
		AutoParameter(const AutoParameter& param);

		/**
		 * @brief Assignment operator.
		 *
		 * @param param The parameter to copy.
		 */
		AutoParameter& operator=(const AutoParameter& param);
	
		virtual ~AutoParameter() {}
	
		AutoParameter* clone() const { return new AutoParameter(* this); }
	
  public:	
	
		/**
		 * @brief Set the value of this parameter.
		 *
		 * This method is redefined so that no constraintException is thrown!
		 * When a Constraint is match, we automatically apply a correct value instead.
		 * This correct value is the nearest limit reached by the value, or a value next to
		 * the limit if the limit is not reachable.
		 *
		 * This allow to perform optimization under constraint whith algorithms that are not
		 * initially built for this.
		 *
		 * @param value the new parameter value.
		 * @throw ConstraintException Never thrown!
		 */
		virtual void setValue(double value) noexcept(false);
	
	public: //Specific method:
		
		/**
		 * @brief Set the message handler for this AutoParameter.
		 *
		 * The message handler keeps all messages that the parameter may send.
		 * The default handler is set to standard output, but you can pass any
		 * ostream object (cerr, ofstream, etc.).
		 *
		 * A NULL pointer disable message output.
		 * 
		 * @param mh The message handler to use.
		 */
		virtual void setMessageHandler(OutputStream* mh) { messageHandler_ = mh; }
		
	public:
	
		static std::string CONSTRAINTS_AUTO;
		static std::string CONSTRAINTS_IGNORE;
		static std::string CONSTRAINTS_KEEP;
};

} //end of namespace bpp.

#endif	//_AUTOPARAMETER_H_

