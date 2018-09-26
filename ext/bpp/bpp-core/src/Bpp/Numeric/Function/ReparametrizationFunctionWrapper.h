//
// File: ReparametrizationFunctionWrapper.h
// Created by: Julien Dutheil
// Created on: Fri Jan  30 09:30 2009
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

#ifndef _REPARAMETRIZATIONFUNCTIONWRAPPER_H_
#define _REPARAMETRIZATIONFUNCTIONWRAPPER_H_

#include "Functions.h"
#include "../AbstractParametrizable.h"
#include "../TransformedParameter.h"

namespace bpp {

/**
 * @brief Function wrapper that remove simple constraints on parameters.
 *
 * This function takes as input another function and reparametrize it when possible.
 * currently, only constraint of the type ]a, b[ where a and b can be +/- inf.
 */
class ReparametrizationFunctionWrapper:
  public virtual Function,
  public AbstractParametrizable
{
   private:
    Function* function_;
    ParameterList functionParameters_;
    
  public:
    /**
     * @brief Build a new reparametrization wrapper for the given function, using all available parameters.
     *
     * @param function The function to reparametrize.
     * @param verbose Print some information.
     */
    ReparametrizationFunctionWrapper(Function* function, bool verbose=true) :
      AbstractParametrizable(function->getNamespace()),
      function_(function),
      functionParameters_(function->getParameters())
    {
      init_(verbose);
    }

    /**
     * @brief Build a new reparametrization wrapper for the given function, using only the specified parameters.
     *
     * @param function The function to reparametrize.
     * @param parameters The list of parameters that will be reparametrized. The intersection with the
     * list of function parameters will be used in the reparametrized function. Any other parameters
     * (in the given list or in the original function) will be ignored.
     * @param verbose Print some information.
     */
    ReparametrizationFunctionWrapper(Function* function, const ParameterList& parameters, bool verbose=true) :
      AbstractParametrizable(function->getNamespace()),
      function_(function),
      functionParameters_(function->getParameters().getCommonParametersWith(parameters))
    {
      init_(verbose);
    }
    
    ReparametrizationFunctionWrapper(const ReparametrizationFunctionWrapper& rfw) :
      AbstractParametrizable(rfw),
      function_(rfw.function_),
      functionParameters_(rfw.functionParameters_) {}
    
    ReparametrizationFunctionWrapper& operator=(const ReparametrizationFunctionWrapper& rfw)
    {
      AbstractParametrizable::operator=(rfw),
      function_ = rfw.function_;
      functionParameters_ = rfw.functionParameters_;
      return *this;
    }

    virtual ~ReparametrizationFunctionWrapper() {}

    ReparametrizationFunctionWrapper* clone() const { return new ReparametrizationFunctionWrapper(*this); }

  private:
    void init_(bool verbose);

  public:
    virtual const Function& getFunction() const { return *function_; } 
    
    virtual Function& getFunction() { return *function_; } 

    void setParameters(const ParameterList& parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
//      parameters.printParameters(std::cout);
      matchParametersValues(parameters);
      //We only set parameters that have been changed:
//      functionParameters_.printParameters(std::cout);
      function_->setParameters(functionParameters_.subList(parameters.getParameterNames()));
    }

    double getValue() const throw (Exception)
    {
      return function_->getValue();
    }

    void fireParameterChanged (const ParameterList &parameters);

};

/**
 * @brief Function wrapper that remove simple constraints on parameters. Also transform first order derivatives.
 *
 * This function takes as input another function and reparametrize it when possible.
 * currently, only constraint of the type ]a, b[ where a and b can be +/- inf.
 */
class ReparametrizationDerivableFirstOrderWrapper:
  public virtual DerivableFirstOrder,
  public ReparametrizationFunctionWrapper
{
    
  public:
    /**
     * @brief Build a new reparametrization wrapper for the given function, using all available parameters.
     *
     * @param function The function to reparametrize.
     * @param verbose Print some information.
     */
    ReparametrizationDerivableFirstOrderWrapper(DerivableFirstOrder* function, bool verbose=true) :
      ReparametrizationFunctionWrapper(function, verbose)
    {}

    /**
     * @brief Build a new reparametrization wrapper for the given function, using only the specified parameters.
     *
     * @param function The function to reparametrize.
     * @param parameters The list of parameters that will be reparametrized. The intersection with the
     * list of function parameters will be used in the reparametrized function. Any other parameters
     * (in the given list or in the original function) will be ignored.
     * @param verbose Print some information.
     */
    ReparametrizationDerivableFirstOrderWrapper(DerivableFirstOrder* function, const ParameterList& parameters, bool verbose=true) :
      ReparametrizationFunctionWrapper(function, parameters, verbose)
    {}
    
    virtual ~ReparametrizationDerivableFirstOrderWrapper() {}

    ReparametrizationDerivableFirstOrderWrapper* clone() const { return new ReparametrizationDerivableFirstOrderWrapper(*this); }

  private:
    void init_(bool verbose);

  public:
    void enableFirstOrderDerivatives(bool yn) { dynamic_cast<DerivableFirstOrder&>(getFunction()).enableFirstOrderDerivatives(yn); }
    
    bool enableFirstOrderDerivatives() const { return dynamic_cast<const DerivableFirstOrder&>(getFunction()).enableFirstOrderDerivatives(); }

    double getFirstOrderDerivative(const std::string& variable) const throw (Exception)
    {
      return dynamic_cast<const DerivableFirstOrder&>(getFunction()).getFirstOrderDerivative(variable)
           * dynamic_cast<const TransformedParameter&>(getParameter(variable)).getFirstOrderDerivative();
    }

};


/**
 * @brief Function wrapper that remove simple constraints on parameters. Also transform first and second order derivatives.
 *
 * This function takes as input another function and reparametrize it when possible.
 * currently, only constraint of the type ]a, b[ where a and b can be +/- inf.
 */
class ReparametrizationDerivableSecondOrderWrapper:
  public virtual DerivableSecondOrder,
  public ReparametrizationDerivableFirstOrderWrapper
{
    
  public:
    /**
     * @brief Build a new reparametrization wrapper for the given function, using all available parameters.
     *
     * @param function The function to reparametrize.
     * @param verbose Print some information.
     */
    ReparametrizationDerivableSecondOrderWrapper(DerivableSecondOrder* function, bool verbose=true) :
      ReparametrizationDerivableFirstOrderWrapper(function, verbose)
    {}

    /**
     * @brief Build a new reparametrization wrapper for the given function, using only the specified parameters.
     *
     * @param function The function to reparametrize.
     * @param parameters The list of parameters that will be reparametrized. The intersection with the
     * list of function parameters will be used in the reparametrized function. Any other parameters
     * (in the given list or in the original function) will be ignored.
     * @param verbose Print some information.
     */
    ReparametrizationDerivableSecondOrderWrapper(DerivableSecondOrder* function, const ParameterList& parameters, bool verbose=true) :
      ReparametrizationDerivableFirstOrderWrapper(function, parameters, verbose)
    {}
    
    virtual ~ReparametrizationDerivableSecondOrderWrapper() {}

    ReparametrizationDerivableSecondOrderWrapper* clone() const { return new ReparametrizationDerivableSecondOrderWrapper(*this); }

  private:
    void init_(bool verbose);

  public:
    void enableSecondOrderDerivatives(bool yn) { dynamic_cast<DerivableSecondOrder&>(getFunction()).enableSecondOrderDerivatives(yn); }
    
    bool enableSecondOrderDerivatives() const { return dynamic_cast<const DerivableSecondOrder&>(getFunction()).enableSecondOrderDerivatives(); }

    double getSecondOrderDerivative(const std::string& variable) const throw (Exception)
    {
      return dynamic_cast<const DerivableSecondOrder&>(getFunction()).getSecondOrderDerivative(variable)
           * pow(dynamic_cast<const TransformedParameter&>(getParameter(variable)).getFirstOrderDerivative(), 2)
           + dynamic_cast<const DerivableSecondOrder&>(getFunction()).getFirstOrderDerivative(variable)
           * dynamic_cast<const TransformedParameter&>(getParameter(variable)).getSecondOrderDerivative();
    }

    double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const throw (Exception)
    {
      return dynamic_cast<const DerivableSecondOrder&>(getFunction()).getSecondOrderDerivative(variable1, variable2)
           * dynamic_cast<const TransformedParameter&>(getParameter(variable1)).getFirstOrderDerivative()
           * dynamic_cast<const TransformedParameter&>(getParameter(variable2)).getFirstOrderDerivative();
    }
 
};

} //end of namespace bpp.

#endif //_REPARAMETRIZATIONFUNCTIONWRAPPER_H_

