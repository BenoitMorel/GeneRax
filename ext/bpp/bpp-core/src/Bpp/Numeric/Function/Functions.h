//
// File: Functions.h
// Created by: Julien Dutheil
// Created on: Sun Nov  9 23:11:00 2003
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

#ifndef _FUNCTIONS_H_
#define _FUNCTIONS_H_

#include "../ParameterList.h"
#include "../Parametrizable.h"
#include "../AbstractParametrizable.h"
#include "../ParameterExceptions.h"

// From Utils:
#include "../../Clonable.h"
#include "../../Exceptions.h"

// From the STL:
#include <cmath>

namespace bpp
{

/**
 * @brief This is the function abstract class.
 *
 * This class provides the interface for function objet
 * and a default implementation of the f() function.
 *
 * The f() function sends the value of the function according to a
 * given set of parameters.
 *
 * However for complexe function like likelihood for instance,
 * computing the function value takes some time, and one do not want
 * to perform several times the computation for an identical set of 
 * parameters.
 * The setParameters() method hence allows to set the parameter value
 * for which the function is to be computed, perform the computation
 * and store the results.
 * The getValue() methods send the result of the computation.
 * One may hence access to the result of the computation by calling the
 * getvalue() method without re-computating the function.
 * The f(parameters) function is a shortcut for
 * @code
 * setParameters(parameters);
 * return getValue();
 * @endcode
 * for convinience.
 *
 * @see Parameter, ParameterList
 */
class Function:
  public virtual Parametrizable
{    
  public:
    Function() {}
    virtual ~Function() {}

  public:

    /**
     * @brief Set the point where the function must be computed.
     *
     * @param parameters The parameter set to pass to the function.
     */
    virtual void setParameters(const ParameterList& parameters)
      noexcept(false) = 0;

    /**
     * @brief Get the value of the function at the current point.
     *
     * @return The value of the function.
     * @throw Exception If no point is specified or if an error occured.
     */
    virtual double getValue() const noexcept(false) = 0;
    
    /**
     * @brief Get the value of the function according to a given set of parameters.
     * 
     * @param parameters The parameter set to pass to the function.
     * @return The value of the function with the given parameter set.
     * @throw Exception If an error occured.
     */
    virtual double f(const ParameterList& parameters) noexcept(false)
    {
      setParameters(parameters);
      return getValue();
    }
};

/**
 * @brief This is the abstract class for first order derivable functions.
 *
 * This class adds the getFirstOrderDerivative() and df() shortcut functions.
 */
class DerivableFirstOrder:
  public virtual Function
{
 public:
    DerivableFirstOrder() {}
    virtual ~DerivableFirstOrder() {}

    DerivableFirstOrder* clone() const = 0;

  public:

    /**
     * @brief Tell if derivatives must be computed.
     *
     * @param yn yes/no
     */
    virtual void enableFirstOrderDerivatives(bool yn) = 0;
    
    /**
     * @brief Tell if derivatives must be computed.
     *
     * @return yes/no
     */
    virtual bool enableFirstOrderDerivatives() const = 0;

    /**
     * @brief Get the derivative of the function at the current point.
     *
     * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{df}{dx} @f$.
     * @return The value of the function.
     * @throw Exception If no point is specified or if an error occured.
     */
    virtual double getFirstOrderDerivative(const std::string& variable) const
      noexcept(false) = 0;
    
    /**
     * @brief Get the value of the first derivative of the function
     * according to a given set of parameters.
     *
     * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{df}{dx} @f$.
     * @param parameters The parameter set to pass to the function.
     * @return The value of the function with the given parameter set.
     * @throw Exception If an error occured.
     */
    virtual double df(const std::string& variable, const ParameterList& parameters)
      noexcept(false)
    {
      setParameters(parameters);
      return getFirstOrderDerivative(variable);
    }
};

/**
 * @brief This is the abstract class for second order derivable functions.
 * 
 * This class adds the getSecondOrderDerivative() and d2f() shortcut functions.
 * Cross derivative functions are also provided.
 */
class DerivableSecondOrder:
  public virtual DerivableFirstOrder
{
  public:
    DerivableSecondOrder() {}
    virtual ~DerivableSecondOrder() {}

    DerivableSecondOrder* clone() const = 0;

  public:

    /**
     * @brief Tell if derivatives must be computed.
     *
     * @param yn yes/no
     */
    virtual void enableSecondOrderDerivatives(bool yn) = 0;
    
    /**
     * @brief Tell if derivatives must be computed.
     *
     * @return yes/no
     */
    virtual bool enableSecondOrderDerivatives() const = 0;

    /**
     * @brief Get the second order derivative of the function at the current point.
     *
     * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x^2} @f$.
     * @return The value of the function.
     * @throw Exception If no point is specified or if an error occured.
     */
    virtual double getSecondOrderDerivative(const std::string& variable) const
      noexcept(false) = 0;
  
    /**
     * @brief Get the value of the second order derivative of the function
     * according to a given set of parameters.
     *
     * @param variable   The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x^2} @f$.
     * @param parameters The parameter set to pass to the function.
     * @return The value of the function with the given parameter set.
     * @throw Exception If an error occured.
     */
    virtual double d2f(const std::string& variable, const ParameterList& parameters)
      noexcept(false)
    {
      setParameters(parameters);
      return getSecondOrderDerivative(variable);
    }    

    /**
     * @brief Get the value of the cross derivative of the function
     * according to a given set of parameters.
     *
     * @param variable1  The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
     * @param variable2  The name of the @f$ y @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
     * @return The value of the function with the given parameter set.
     * @throw Exception If an error occured.
     */
    virtual double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const
      noexcept(false) = 0;
    
    /**
     * @brief Get the value of the cross derivative of the function
     * according to a given set of parameters.
     *
     * @param variable1  The name of the @f$ x @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
     * @param variable2  The name of the @f$ y @f$ variable in @f$ \frac{\partial^2 f}{\partial x \partial y} @f$.
     * @param parameters The parameter set to pass to the function.
     * @return The value of the function with the given parameter set.
     * @throw Exception If an error occured.
     */
    virtual double d2f(const std::string& variable1, const std::string& variable2, const ParameterList& parameters)
      noexcept(false)
    {
      setParameters(parameters);
      return getSecondOrderDerivative(variable1, variable2);
    }
};

/**
 * @brief General class that wraps a function into another one.
 * This class is meant to be derivated and just provides a general framework.
 */
class FunctionWrapper:
  public virtual Function
{
  protected:
    Function* function_;

  public:
    FunctionWrapper(Function* function) : function_(function) {}
    FunctionWrapper(const FunctionWrapper& fw) : function_(fw.function_) {}
    FunctionWrapper& operator=(const FunctionWrapper& fw)
    {
      function_ = fw.function_;
      return *this;
    }

  public:
    bool hasParameter(const std::string& name) const
    {
      return function_->hasParameter(name);
    }

    void setParameters(const ParameterList & parameters)
      noexcept(false)
    {
      function_->setParameters(parameters);
    }

    const ParameterList& getParameters() const noexcept(false)
    {
      return function_->getParameters();  
    }

    const Parameter& getParameter(const std::string & name) const
      noexcept(false)
    {
      return function_->getParameter(name);
    }

    double getValue() const noexcept(false)
    {
      return function_->getValue();
    }
    
    double f(const ParameterList& parameters) noexcept(false)
    {
      return function_->f(parameters);
    }
    
    double getParameterValue(const std::string& name) const noexcept(false)
    {
      return function_->getParameterValue(name);
    }
      
    void setAllParametersValues(const ParameterList & parameters)
      noexcept(false)
    {
      function_->setAllParametersValues(parameters);
    }
    
    void setParameterValue(const std::string& name, double value)
      noexcept(false)
    {
      function_->setParameterValue(name, value);
    }
    
    void setParametersValues(const ParameterList& parameters)
      noexcept(false)
    {
      function_->setParametersValues(parameters);
    }
    
    bool matchParametersValues(const ParameterList& parameters)
      noexcept(false)
    {
      return function_->matchParametersValues(parameters);
    }

    size_t getNumberOfParameters() const
    {
      return function_->getNumberOfParameters();
    }

    void setNamespace(const std::string& prefix)
    {
      function_->setNamespace(prefix);
    }

    std::string getNamespace() const
    {
      return function_->getNamespace();
    }

    std::string getParameterNameWithoutNamespace(const std::string& name) const
    {
      return function_->getParameterNameWithoutNamespace(name);
    }

};



/**
 * @brief General class that wraps a function into another one.
 * This class is meant to be derivated and just provides a general framework.
 */
class DerivableFirstOrderWrapper:
  public FunctionWrapper,
  public virtual DerivableFirstOrder
{
  public:
    DerivableFirstOrderWrapper(DerivableFirstOrder* function) : FunctionWrapper(function) {}

  public:
    void enableFirstOrderDerivatives(bool yn) {
      dynamic_cast<DerivableFirstOrder*>(function_)->enableFirstOrderDerivatives(yn);
    }
    
    bool enableFirstOrderDerivatives() const {
      return dynamic_cast<DerivableFirstOrder*>(function_)->enableFirstOrderDerivatives();
    }

    double getFirstOrderDerivative(const std::string& variable) const
      noexcept(false) {
      return dynamic_cast<DerivableFirstOrder*>(function_)->getFirstOrderDerivative(variable);
    }

};



/**
 * @brief General class that wraps a function into another one.
 * This class is meant to be derivated and just provides a general framework.
 */
class DerivableSecondOrderWrapper:
  public DerivableFirstOrderWrapper,
  public virtual DerivableSecondOrder
{
  public:
    DerivableSecondOrderWrapper(DerivableSecondOrder* function) : DerivableFirstOrderWrapper(function) {}

  public:
    void enableSecondOrderDerivatives(bool yn) {
      dynamic_cast<DerivableSecondOrder*>(function_)->enableSecondOrderDerivatives(yn);
    }
    
    bool enableSecondOrderDerivatives() const {
      return dynamic_cast<DerivableSecondOrder*>(function_)->enableSecondOrderDerivatives();
    }

    double getSecondOrderDerivative(const std::string& variable) const
      noexcept(false) {
      return dynamic_cast<DerivableSecondOrder*>(function_)->getSecondOrderDerivative(variable);
    }

    double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const
      noexcept(false) {
      return dynamic_cast<DerivableSecondOrder*>(function_)->getSecondOrderDerivative(variable1, variable2);
    }

};



/**
 * @brief Wrapper class for optimization under constraints.
 *
 * Catch any ConstraintException thrown and send +inf.
 */
class InfinityFunctionWrapper:
  public FunctionWrapper
{
  protected:
    mutable bool constraintMatch_;
    
  public:
    InfinityFunctionWrapper(Function* function) :
      FunctionWrapper(function),
      constraintMatch_(false) {}
    virtual ~InfinityFunctionWrapper() {}

    InfinityFunctionWrapper* clone() const { return new InfinityFunctionWrapper(*this); }

  public:

    void setParameters(const ParameterList& parameters)
    noexcept(false)
    {
      try
      {
        function_->setParameters(parameters);
        constraintMatch_ = false;
      }
      catch(ConstraintException& ce)
      {
        constraintMatch_ = true;
      }
    }

    double getValue() const
    noexcept(false)
    {
      return constraintMatch_ ? -log(0.) :  function_->getValue();
    }
    
    double f(const ParameterList& parameters) noexcept(false)
    {
      setParameters(parameters);
      return getValue();
    }
          
    void setAllParametersValues(const ParameterList & parameters)
    noexcept(false)
    {
      try
      {
        function_->setAllParametersValues(parameters);
        constraintMatch_ = false;
      }
      catch(ConstraintException& ce)
      {
        constraintMatch_ = true;
      }
    }
    
    void setParameterValue(const std::string& name, double value)
    noexcept(false)
    {
      try
      {
        function_->setParameterValue(name, value);
        constraintMatch_ = false;
      }
      catch(ConstraintException& ce)
      {
        constraintMatch_ = true;
      }
    }
    
    void setParametersValues(const ParameterList& parameters)
    noexcept(false)
    {
      try
      {
        function_->setParametersValues(parameters);
        constraintMatch_ = false;
      }
      catch(ConstraintException& ce)
      {
        constraintMatch_ = true;
      }
    }
    
    bool matchParametersValues(const ParameterList& parameters)
    noexcept(false)
    {
      try
      {
        bool test = function_->matchParametersValues(parameters);
        constraintMatch_ = false;
        return test;
      }
      catch (ConstraintException& ce)
      {
        constraintMatch_ = true;
        return false;
      }
    }

};

/**
 * @brief Wrapper class for optimization under constraints.
 *
 * Catch any ConstraintException thrown and send +inf.
 */
class InfinityDerivableFirstOrderWrapper :
  public virtual InfinityFunctionWrapper
{
  public:
    InfinityDerivableFirstOrderWrapper(DerivableFirstOrder* function) : InfinityFunctionWrapper(function) {}
    virtual ~InfinityDerivableFirstOrderWrapper() {}
    
    InfinityDerivableFirstOrderWrapper* clone() const { return new InfinityDerivableFirstOrderWrapper(*this); }

  public:
    
    double getFirstOrderDerivative(const std::string& variable) const
    noexcept(false)
    {
      return constraintMatch_ ? -log(0.) :  (dynamic_cast<DerivableFirstOrder *>(function_)->getFirstOrderDerivative(variable));    
    }
    
    double df(const std::string& variable, const ParameterList& parameters)
    noexcept(false)
    {
      setParameters(parameters);
      return getFirstOrderDerivative(variable);
    }
};

/**
 * @brief Wrapper class for optimization under constraints.
 *
 * Catch any ConstraintException thrown and send +inf.
 */
class InfinityDerivableSecondOrderWrapper :
  public virtual InfinityDerivableFirstOrderWrapper
{
  public:
    InfinityDerivableSecondOrderWrapper(DerivableFirstOrder* function):
      InfinityFunctionWrapper(function),
      InfinityDerivableFirstOrderWrapper(function) {}
    virtual ~InfinityDerivableSecondOrderWrapper() {}

    InfinityDerivableSecondOrderWrapper* clone() const { return new InfinityDerivableSecondOrderWrapper(*this); }

  public:

    double getSecondOrderDerivative(const std::string& variable) const
    noexcept(false)
    {
      return constraintMatch_ ? -log(0.) :  (dynamic_cast<DerivableSecondOrder *>(function_)->getSecondOrderDerivative(variable));          
    }
  
    double d2f(const std::string & variable, const ParameterList& parameters)
    noexcept(false)
    {
      setParameters(parameters);
      return getSecondOrderDerivative(variable);
    }    

    double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const
    noexcept(false)
    {
      return constraintMatch_ ? -log(0.) :  (dynamic_cast<DerivableSecondOrder *>(function_)->getSecondOrderDerivative(variable1, variable2));      
    }
    
    double d2f(const std::string & variable1, const std::string& variable2, const ParameterList& parameters)
    noexcept(false)
    {
      setParameters(parameters);
      return getSecondOrderDerivative(variable1, variable2);
    }
};


/**
 * @brief A simple funciton with two parameters, mostly for testing and debugging :)
 *
 * @author Julien Dutheil.
 */
class TestFunction :
  public virtual Function,
  public AbstractParametrizable
{
  public:
    TestFunction(double x = 0, double y = 0) :
      AbstractParametrizable("")
    {
      addParameter_(new Parameter("x", x));
      addParameter_(new Parameter("y", y));
    }

    Clonable* clone() const { return new TestFunction(*this); }

    void setParameters(const ParameterList& parameters) noexcept(false)
    {
      matchParametersValues(parameters);
    }

    double getValue() const noexcept(false)
    {
      double x = getParameter("x").getValue();
      double y = getParameter("y").getValue();
      return (x*x + y*y);
    }

    void fireParameterChanged(const ParameterList& parameters) {}
};

} //end of namespace bpp.

#endif  //_FUNCTIONS_H_

