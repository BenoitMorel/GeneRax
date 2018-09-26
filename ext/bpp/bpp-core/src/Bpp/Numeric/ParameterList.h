//
// File: ParameterList.h
// Created by: Julien Dutheil
// Created on: Wed Oct 15 18:17:29 2003
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

#ifndef _PARAMETERLIST_H_
#define _PARAMETERLIST_H_

#include "Parameter.h"
#include "../Clonable.h"
#include "../Io/OutputStream.h"

// From STL:
#include <vector>
#include <string>
#include <iostream>

namespace bpp
{
/**
 * @brief The parameter list object.
 *
 * @author Julien Dutheil, Laurent Gueguen
 * This is a vector of Parameter with a few additional methods, mainly for giving
 * name access.
 */
class ParameterList :
  public Clonable
{
private:
  std::vector<std::shared_ptr<Parameter> > parameters_;

public:
  /**
   * @brief Build a new ParameterList object.
   */
  ParameterList() : parameters_() {}

  /**
   * @brief Copy constructor
   *
   * All parameters in the list will be cloned.
   */
  ParameterList(const ParameterList& pl);

  ParameterList& operator=(const ParameterList& pl);

  ParameterList* clone() const { return new ParameterList(*this); }

  virtual ~ParameterList();

public:
  /**
   * @return The number of parameters in the list.
   */
  size_t size() const { return parameters_.size(); }

  /**
   * @return The parameter at a given position.
   * @warning No check is performed on the validity of the index given as input!
   */
  virtual const Parameter& operator[](size_t i) const { return *parameters_[i]; }
  virtual Parameter& operator[](size_t i) { return *parameters_[i]; }

  /**
   * @return The shared_ptr parameter at a given position.
   * @warning No check is performed on the validity of the index given as input!
   */
  virtual const std::shared_ptr<Parameter>& getSharedParameter(size_t i) const { return parameters_[i]; }
  virtual std::shared_ptr<Parameter>& getSharedParameter(size_t i) { return parameters_[i]; }

  /**
   * @brief Get the parameter with name <i>name</i>.
   *
   * @param name The name of the parameter to look for.
   * @return A const reference toward the parameter with name <i>name</i>.
   * @throw ParameterNotFoundException If no parameter with the given name is found.
   */
  virtual const Parameter& getParameter(const std::string& name) const noexcept(false);

  /**
   * @brief Get the parameter with name <i>name</i> as a shared pointer.
   *
   * @param name The name of the parameter to look for.
   * @return A const shared parameter toward the parameter with name <i>name</i>.
   * @throw ParameterNotFoundException If no parameter with the given name is found.
   */
  virtual const std::shared_ptr<Parameter>& getSharedParameter(const std::string& name) const noexcept(false);

  /**
   * @brief Get the value of the parameter with name <i>name</i>.
   *
   * @param name The name of the parameter to look for.
   * @return A value of the parameter with name <i>name</i>.
   * @throw ParameterNotFoundException If no parameter with the given name is found.
   */

  virtual double getParameterValue(const std::string& name) const noexcept(false);

  /**
   * @brief Get the parameter with name <i>name</i>.
   *
   * @param name The name of the parameter to look for.
   * @return A reference toward the parameter with name <i>name</i>.
   * @throw ParameterNotFoundException If no parameter with the given name is found.
   */
  virtual Parameter& getParameter(const std::string& name) noexcept(false);

  /**
   * @brief Get the parameter with name <i>name</i> as a shared pointer.
   *
   * @param name The name of the parameter to look for.
   * @return A shared parameter toward the parameter with name <i>name</i>.
   * @throw ParameterNotFoundException If no parameter with the given name is found.
   */
  virtual std::shared_ptr<Parameter>& getSharedParameter(const std::string& name) noexcept(false);

  /**
   * @brief Get given parameters as a sublist.
   *
   * @param names Name of the parameters to be included in the list.
   * @return A list with all parameters specified.
   * @throw ParameterNotFoundException If at least one name does not correspond to a parameter in the list.
   */
  virtual ParameterList subList(const std::vector<std::string>& names) const noexcept(false);

  /**
   * @brief Get given parameter as a sublist.
   *
   * @param name Name of the parameter to be included in the list.
   * @return A list with the parameter specified.
   * @throw ParameterNotFoundException If no parameter with the given name is found.
   */
  virtual ParameterList subList(const std::string& name) const noexcept(false);

  /**
   * @brief Get given parameters as a sublist.
   *
   * @param parameters Positions of the parameters to be included in the list.
   * @return A list with all parameters specified.
   */
  virtual ParameterList subList(const std::vector<size_t>& parameters) const;

  /**
   * @brief Get given parameter as a sublist.
   *
   * @param parameter Position of the parameters to be included in the list.
   * @return A list with the parameter specified.
   */
  virtual ParameterList subList(size_t parameter) const;

  /**
   * @brief Get the sublist containing all common parameter between this list and pl.
   *
   * @param params The list to compare to.
   * @return A list with all common parameters.
   */
  virtual ParameterList getCommonParametersWith(const ParameterList& params) const;

  /**
   * @brief Get all parameter names in the list.
   *
   * @return A vector with all names in the same order as the parameters in the list.
   */
  virtual std::vector<std::string> getParameterNames() const;

  /**
   * @brief Get all parameter names matching with the given name. Up
   * to now, only "*" jokers are available.
   *
   * @param pattern a pattern of name
   * @return A vector of matching names
   */
  
  virtual std::vector<std::string> getMatchingParameterNames(const std::string& pattern) const;

  /**
   * @brief Add a new parameter at the end of the list.
   *
   * @param param The parameter to add to the list.
   */
  virtual void addParameter(const Parameter& param) noexcept(false);

  /**
   * @brief Add a new parameter at the end of the list.
   *
   * This function is more efficient than its reference counterpart,
   * as it avoid an object copy. The ParameterList will own the pointer
   * after addition, if it is successful.
   *
   * @param param A ppointer toward the parameter to add to the list.
   */

  virtual void addParameter(Parameter* param) noexcept(false);

  /**
   * @brief Share a parameter at the end of the list.
   *
   * @param param The shared_ptr parameter to add to the list.
   */
  
  virtual void shareParameter(const std::shared_ptr<Parameter>& param) noexcept(false);


  /**
   * @brief Change given parameter.
   *
   * @param index The position of the parameter to alter.
   * @param param The parameter to add to the list.
   * @throw IndexOutOfBoundsException if the index is not valid.
   */
  
  virtual void setParameter(size_t index, const Parameter& param) noexcept(false);

//  virtual void setParameter(size_t index, Parameter* param) throw (IndexOutOfBoundsException);

  /**
   * @brief Add new parameters at the end of the list.
   *
   * @param params The parameter list containing the new parameters to
   * add to the list.
   */
  virtual void addParameters(const ParameterList& params) noexcept(false);

  /**
   * @brief Share parameters with a given list. They are added the end of this list.
   *
   * @param params The parameter list containing the parameters to
   * share with to the list.
   */
  
  virtual void shareParameters(const ParameterList& params) noexcept(false);

  /**
   * @brief Add parameters to the list. If the parameter already
   * exists, only the value is updated, otherwise the new parameter is
   * added at the end of the list.
   *
   * @param params The parameter list containing the new paramters to add to the list.
   */
  virtual void includeParameters(const ParameterList& params);

  /**
   * @brief Set the value of parameter with name <i>name</i> to be equal to <i>value</i>.
   *
   * @param name the name of the parameter to set.
   * @param value The value of the parameter.
   * @throw ParameterNotFoundException If no parameter with the given name is found in the list.
   * @throw ConstraintException If the value is incorrect.
   */
  virtual void setParameterValue(const std::string& name, double value)
  noexcept(false);

  /**
   * @brief Set the parameters to be equals to <i>params</i>.
   *
   * The list must contain exactly the same parameters (ie same names)
   * than the parameters available.
   *
   * @param params A list with all parameters.
   * @see setParameters(), matchParameters();
   * @throw ParameterNotFoundException If at least one name does not correspond to a parameter in the list.
   * @throw ConstraintException If one value is incorrect (and the two parameter list do not have the same constraints).
   */
  virtual void setAllParametersValues(const ParameterList& params)
  noexcept(false);

  /**
   * @brief Update the parameters from the ones in <i>params</i>
   * that have matching names.
   *
   * @param params A list containing all parameters to update.
   * @see setAllParameters(), matchParameters()
   * @throw ConstraintException If one value is incorrect (and the two parameter list do not have the same constraints).
   */
  virtual void setParametersValues(const ParameterList& params);

  /**
   * @brief Returns true if the Parameter of the given name exists.
   *
   * @name A string name
   */
  virtual bool hasParameter(const std::string& name) const;

  /**
   * @brief Tests the parameters from <i>params</i>.
   *
   * Only common parameters with <i>params</i> are compared.
   *
   * @param params A list of parameters.
   * @return true iff a least one parameter value is different.
   */
  virtual bool testParametersValues(const ParameterList& params) const;

  /**
   * @brief Update the parameters from <i>params</i>.
   *
   * Only common parameters with <i>params</i> will be updated.
   *
   * @param params A list of parameters.
   * @param updatedParameters An optional pointer toward a vector which will
   * store the indices of parameters for which a value has changed.
   * Indices are relative on the input parameter list "params".
   * @return true iff a least one parameter value has been changed.
   * @see setParameters(), setAllParameters()
   */

  virtual bool matchParametersValues(const ParameterList& params, std::vector<size_t>* updatedParameters = 0)
  noexcept(false);

  /**
   * @brief Set the parameters to be equals to <i>params</i>.
   *
   * The list must contain exactly the same parameters (ie same names)
   * than the parameters available.
   *
   * @param params A list with all parameters.
   * @see setParameters(), matchParameters();
   */
  virtual void setAllParameters(const ParameterList& params)
  noexcept(false);

  /**
   * @brief Update the parameters from <i>params</i>.
   *
   * <i>params</i> must be a subset of all parameters available.
   *
   * @param params A list containing all parameters to update.
   * @see setAllParameters(), matchParameters()
   */
  virtual void setParameters(const ParameterList& params)
  noexcept(false);

  /**
   * @brief Update the parameters from <i>params</i>.
   *
   * Only common parameters with <i>params</i> will be updated.
   *
   * @param params A list of parameters.
   * @see setParameters(), setAllParameters()
   */
  virtual void matchParameters(const ParameterList& params);

  /**
   * @brief Delete a parameter from the list.
   *
   * @param name The name of the parameter to delete from the list.
   */
  virtual void deleteParameter(const std::string& name) noexcept(false);

  /**
   * @brief Delete several parameters from the list.
   *
   * @param names The names of the parameters to delete from the list.
   * @param mustExist If true, an exception is thrown if a name does
   *                    not match an extant parameter.
   */
  
  virtual void deleteParameters(const std::vector<std::string>& names,
                                bool mustExist = true)
  noexcept(false);

  /**
   * @brief Delete a parameter from the list.
   *
   * @param index The position of the parameter to delete in the list.
   */
  virtual void deleteParameter(size_t index) noexcept(false);

  /**
   * @brief Delete several parameters from the list.
   *
   * @param indices The positions of the parameters to delete in the list.
   * Duplicated positions will be considered only one time.
   */
  virtual void deleteParameters(const std::vector<size_t>& indices)
  noexcept(false);

  /**
   * @brief Get the position of a given parameter according to its name.
   *
   * @param name The name of the parameter to look for.
   * @return The position of the parameter if found. If several parameters exist with the given name,
   * the position of the first one is returned.
   * @throw ParameterNotFoundException If no parameter with the given name is found.
   */
  virtual size_t whichParameterHasName(const std::string& name) const
  noexcept(false);

  /**
   * @brief Print all parameters.
   */
  virtual void printParameters(OutputStream& out) const;

  virtual void printParameters(std::ostream& out) const
  {
    StlOutputStreamWrapper os(&out);
    printParameters(os);
  }

  /**
   * @brief Reset the list: delete all parameters.
   */
  virtual void reset();
};
} // end of namespace bpp.

#endif  // _PARAMETERLIST_H_

