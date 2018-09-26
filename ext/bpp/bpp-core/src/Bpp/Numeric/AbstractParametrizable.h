//
// File: AbstractParametrizable.h
// Created by: Julien Dutheil
// Created on: Sun Mar 29 09:10 2009
// Created from file Parametrizable.h
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

#ifndef _ABSTRACTPARAMETRIZABLE_H_
#define _ABSTRACTPARAMETRIZABLE_H_

#include "Parametrizable.h"

//From the STL:
#include <map>

namespace bpp
{

/**
 * @brief A partial implementation of the Parametrizable interface.
 *
 * Parameters are stored in a protected ParameterList object.
 *
 * The abstract fireParameterChanged() method is provided so that the derived class
 * know when a parameter has changed, and can be updated.
 * All methods call the corresponding method in ParameterList and then call the
 * fireParameterChanged() method.
 */
class AbstractParametrizable:
  public virtual Parametrizable
{
  private:
    ParameterList parameters_;
    std::string prefix_;

  public:
    AbstractParametrizable(const std::string& prefix) : parameters_(), prefix_(prefix) {}

    virtual ~AbstractParametrizable() {}

  public:
    bool hasParameter(const std::string& name) const { return parameters_.hasParameter(prefix_ + name); }

    const ParameterList& getParameters() const { return parameters_; }
    
    const Parameter& getParameter(const std::string& name) const noexcept(false)
    {
      return parameters_.getParameter(prefix_ + name);
    }

  const std::shared_ptr<Parameter>& getSharedParameter(const std::string& name) const
    noexcept(false)
  {
    return parameters_.getSharedParameter(prefix_ + name);
  }

    double getParameterValue(const std::string& name) const
      noexcept(false)
    { 
      return getParameter(name).getValue();
    }

    void setAllParametersValues(const ParameterList & parameters)
      noexcept(false)
    {
      parameters_.setAllParametersValues(parameters);
      fireParameterChanged(parameters);
    }

    void setParameterValue(const std::string& name, double value)
      noexcept(false)
    {
      parameters_.setParameterValue(prefix_ + name, value);
      fireParameterChanged(parameters_.subList(prefix_ + name));
    }

    void setParametersValues(const ParameterList& parameters)
      noexcept(false)
    { 
      parameters_.setParametersValues(parameters);
      fireParameterChanged(parameters);
    }

    bool matchParametersValues(const ParameterList& parameters)
      noexcept(false)
    {
      std::unique_ptr< std::vector<size_t> >updatedParameters(new std::vector<size_t>());
      bool test = parameters_.matchParametersValues(parameters, updatedParameters.get());
      if (test) 
        fireParameterChanged(parameters.subList(*updatedParameters));
      return test;
    }

    size_t getNumberOfParameters() const { return parameters_.size(); }
     
    void setNamespace(const std::string& prefix);
    
    std::string getNamespace() const { return prefix_; }
    
    std::string getParameterNameWithoutNamespace(const std::string& name) const;

    /**
     * @brief Notify the class when one or several parameters have changed.
     *
     * @param parameters A ParameterList object with parameters that changed.
     */
  
  virtual void fireParameterChanged(const ParameterList& parameters) {};
  
  

  protected:
  
  virtual void addParameter_(Parameter* parameter)
  {
    if (parameter)
      parameters_.addParameter(parameter);
  }
  
  virtual void addParameters_(const ParameterList& parameters)
  {
    parameters_.addParameters(parameters);
  }

  void shareParameter_(const std::shared_ptr<Parameter>& parameter)
  {
    parameters_.shareParameter(parameter);
  }

  void shareParameters_(const ParameterList& parameters)
  {
    parameters_.shareParameters(parameters);
  }

  virtual void includeParameters_(const ParameterList& parameters)
  {
    parameters_.includeParameters(parameters);
  }

  virtual void deleteParameter_(size_t index) noexcept(false)
  {
    if (index >= parameters_.size())
      throw IndexOutOfBoundsException("AbstractParametrizable::deleteParameter_.", index, 0, parameters_.size() - 1);
    parameters_.deleteParameter(index);
  }

  virtual void deleteParameter_(std::string& name)
  {
    parameters_.deleteParameter(name);
  }

  virtual void deleteParameters_(const std::vector<std::string>& names)
  {
     parameters_.deleteParameters(names);
  }

  void resetParameters_()
  {
    parameters_.reset();
  }

    /**
     * @param name The name of the parameter.
     * @return A reference toward the corresponding parameter.
     * @throw ParameterNotFoundException If no parameter with that name is found in the list.
     */
    Parameter& getParameter_(const std::string& name) noexcept(false)
    {
      return parameters_.getParameter(prefix_ + name);
    }
  
    /**
     * @param name The name of the parameter, including its namespace.
     * @return A reference toward the corresponding parameter.
     * @throw ParameterNotFoundException If no parameter with that name is found in the list.
     */
    Parameter& getParameterWithNamespace_(const std::string& name)
      noexcept(false)
    {
      return getParameter_(name);
    }
    /**
     * @param name The name of the parameter, including its namespace.
     * @return A reference toward the corresponding parameter.
     * @throw ParameterNotFoundException If no parameter with that name is found in the list.
     */
    const Parameter& getParameterWithNamespace_(const std::string& name) const
      noexcept(false)
    {
      return getParameter(name);
    }

    Parameter& getParameter_(size_t index) noexcept(false)
    {
      if (index >= parameters_.size())
        throw IndexOutOfBoundsException("AbstractParametrizable::getParameter_.", index, 0, parameters_.size() - 1);
      return parameters_[index];
    }

  const Parameter& getParameter_(size_t index) const noexcept(false)
    {
      if(index >= parameters_.size())
        throw IndexOutOfBoundsException("AbstractParametrizable::getParameter_.", index, 0, parameters_.size() - 1);
      return parameters_[index];
    }


  ParameterList& getParameters_() { return parameters_; }

  /**
   * @return The shared_ptr parameter at a given position.
   * @warning No check is performed on the validity of the index given as input!
   */ 
  const std::shared_ptr<Parameter>& getSharedParameter(size_t i) const
  {
    return parameters_.getSharedParameter(i);
  }
  
  std::shared_ptr<Parameter>& getSharedParameter(size_t i)
  {
    return parameters_.getSharedParameter(i);
  }

};

} //end of namespace bpp.

#endif //_ABSTRACTPARAMETRIZABLE_H_

