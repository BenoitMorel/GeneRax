//
// File: Parametrizable.h
// Created by: Julien Dutheil
// Created on: Sun Oct 19 23:06:42 2003
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

#ifndef _PARAMETRIZABLE_H_
#define _PARAMETRIZABLE_H_

// From Utils:
#include "../Clonable.h"

// From the STL:
#include <string>

#include "ParameterList.h"

namespace bpp
{

/**
 * @brief This is the interface for all objects that imply parameters.
 *
 * @see Parameter, ParameterList
 */
class Parametrizable:
  public virtual Clonable
{
	public:
		Parametrizable() {}
		virtual ~Parametrizable() {}

	public:
    /**
     * @brief Tell if there is a parameter with specified name.
     *
     * @param name The name of the parameter to look for.
     * @return y/n.
     */
    virtual bool hasParameter(const std::string& name) const = 0;
	
		/**
		 * @brief Get all parameters available.
     *
     * @see getIndependentParameters if some parameters are aliased.
		 * @return A list with all parameters available.
		 */
		virtual const ParameterList& getParameters() const = 0;

    /**
     * @brief Get the parameter with specified name.
     *
     * @param name The name of the parameter to look for.
     * @return The parameter with given name.
     * @throw ParameterNotFoundException if no parameter with this name is found.
     */
    virtual const Parameter& getParameter(const std::string& name) const
    noexcept(false) = 0;
	
		/**
		 * @brief Get the value for parameter of name 'name'.
		 *
		 * @param name The name of the parameter.
		 * @return the value of parameter <i>name</i>.
		 */
		virtual double getParameterValue(const std::string& name) const
    noexcept(false) = 0;

		/**
		 * @brief Set the parameters values to be equals to those of <i>parameters</i>.
		 *
		 * The list must contain exactly the same parameters (ie same names)
		 * than the parameters available.
		 *
		 * @param parameters A list with all parameters.
		 * @throw ParameterNotFoundException If a some parameter in the list is not in <i>params</i>.
		 * @throw ConstraintException If a value in <i>parameters</i> does not match the constraint in the
		 * corresponding parameter in the list.
		 */
		virtual void setAllParametersValues(const ParameterList& parameters)
    noexcept(false) = 0;

		/**
		 * @brief Set the value of parameter with name <i>name</i> to be equal to <i>value</i>.
		 *
		 * @param name the name of the parameter to set.
		 * @param value The value of the parameter.
		 * @throw ParameterNotFoundException If no parameter in the list has the name <i>name</i>.
		 * @throw ConstraintException If <i>value</i> does not match the constraint associated to
		 * parameter <i>name</i>.
		 */
		virtual void setParameterValue(const std::string& name, double value)
    noexcept(false) = 0;

		/**
		 * @brief Update the parameters from <i>parameters</i>.
		 *
		 * <i>parameters</i> must be a subset of all parameters available.
		 *
		 * @param parameters A list containing all parameters to update.
		 * @throw ParameterNotFoundException If a some parameter in <i>params</i> is not in the list.
		 * @throw ConstraintException If a value in <i>parameters</i> does not match the constraint in the
		 * corresponding parameter in the list.
		 */
		virtual void setParametersValues(const ParameterList& parameters)
    noexcept(false) = 0;

		/**
		 * @brief Update the parameters from <i>parameters</i>.
		 *
		 * Only common parameters with <i>parameters</i> will be updated.
		 *
		 * @param parameters A list of parameters.
     * @return True if at least one parameter value has been changed.
		 * @throw ConstraintException If a value in <i>parameters</i> does not match the constraint in the
		 * corresponding parameter in the list.
		 */
		virtual bool matchParametersValues(const ParameterList& parameters)
    noexcept(false) = 0;

    /**
     * @brief Get the number of parameters.
     *
     * @see getNumberOfIndependentParameters If some parameters are aliased.
     * @return The number of parameters.
     */
    virtual size_t getNumberOfParameters() const = 0;

    /**
     * @brief Set the namespace for the parameter names.
     *
     * @param prefix The 'namespace', that is a prefix to add to all parameter names.
     * If parameter names are already prefixed, the new prefix will be used instead.
     */
    virtual void setNamespace(const std::string& prefix) = 0;

    /**
     * @return The current namespace used. This is an empty string if no namespace is currently defined.
     */
    virtual std::string getNamespace() const = 0;

    /**
     * @brief Resolves a parameter name according to the current namespace.
     *
     * @return The parameter name without the namespace prefix, if any.
     */
    virtual std::string getParameterNameWithoutNamespace(const std::string& name) const = 0;

};

/**
 * @brief A low-level implementation of the Parametrizable interface with void functions.
 *
 * @see Parameter, ParameterList, Parametrizable
 */
class ParametrizableAdapter:
  public virtual Parametrizable
{
  protected:
    ParameterList parameters_;
    Parameter parameter_;

	public:
		ParametrizableAdapter(): parameters_(), parameter_() {}
		virtual ~ParametrizableAdapter() {}

	public:

		/**
		 * @name The Parametrizable interface.
		 *
		 * @{
		 */
    bool hasParameter(const std::string & name) const { return parameters_.hasParameter(name); }
		const ParameterList & getParameters() const { return parameters_; }
    const Parameter & getParameter(const std::string & name) const
      noexcept(false) { return parameter_; }
		double getParameterValue(const std::string & name) const
      noexcept(false) { return 0; };
		void setAllParametersValues(const ParameterList & parameters)
      noexcept(false) {}
		void setParameterValue(const std::string & name, double value)
      noexcept(false) {}
		void setParametersValues(const ParameterList & parameters)
      noexcept(false) {}
		bool matchParametersValues(const ParameterList & parameters)
      noexcept(false) { return false ;}
    size_t getNumberOfParameters() const{ return 0; }
    void setNamespace(const std::string& prefix) {}
    std::string getNamespace() const { return ""; }
    std::string getParameterNameWithoutNamespace(const std::string& name) const { return name; }
		/** @} */

};

} //end of namespace bpp.

#endif	//_PARAMETRIZABLE_H_

