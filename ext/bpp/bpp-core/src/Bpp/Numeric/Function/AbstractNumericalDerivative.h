//
// File: AbstractNumericalDerivative.h
// Created by: Julien Dutheil
// Created on: Thu Aug 17 15:00 2006
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

#ifndef _ABSTRACTNUMERICALDERIVATIVE_H_
#define _ABSTRACTNUMERICALDERIVATIVE_H_

#include "Functions.h"
#include "../Matrix/Matrix.h"

//From the STL:
#include <map>
#include <vector>
#include <string>

namespace bpp
{

/**
 * @brief Numerical derivative function wrapper, partial implementation.
 *
 * This class provides a wrapper for Function object, implementing the DerivableSecondOrder interface
 * (Although no cross derivative is implemented for now).
 * Derivations of this class can be used as full DerivableSecondOrder objects, with derivative functions.
 * 
 * Three kinds of constructors are provided: one with a Function object, another with a DerivableFirstOrder object, and one with a DerivableSecondOrder object.
 * In the first case, all derivatives will be computed numerically.
 * In the second case, first order derivative will be computed numerically only if no appropriate analytical derivative is available, second order derivative will always be computed numerically.
 * In the last case, first and second order derivative will be computed numerically only if no appropriate analytical derivative is available.
 */
class AbstractNumericalDerivative:
  public DerivableSecondOrder,
  public FunctionWrapper
{
  protected:
    DerivableFirstOrder *function1_;
    DerivableSecondOrder *function2_;
    double h_;
    std::vector<std::string> variables_;
    mutable std::map<std::string, size_t> index_; //Store positions in array corresponding to variable names.
    std::vector<double> der1_;
    std::vector<double> der2_;
    RowMatrix<double> crossDer2_;
    bool computeD1_, computeD2_, computeCrossD2_;
    
  public:
    AbstractNumericalDerivative(Function* function):
      FunctionWrapper(function), function1_(0), function2_(0),
      h_(0.0001), variables_(), index_(), der1_(), der2_(), crossDer2_(),
      computeD1_(true), computeD2_(true), computeCrossD2_(false) {}

    AbstractNumericalDerivative(DerivableFirstOrder* function):
      FunctionWrapper(function), function1_(function), function2_(0),
      h_(0.0001), variables_(), index_(), der1_(), der2_(), crossDer2_(),
      computeD1_(true), computeD2_(true), computeCrossD2_(false) {}

    AbstractNumericalDerivative(DerivableSecondOrder* function):
      FunctionWrapper(function), function1_(function), function2_(function),
      h_(0.0001), variables_(), index_(), der1_(), der2_(), crossDer2_(),
      computeD1_(true), computeD2_(true), computeCrossD2_(false) {}

    AbstractNumericalDerivative(const AbstractNumericalDerivative& ad):
      FunctionWrapper(ad), function1_(ad.function1_), function2_(ad.function2_),
      h_(ad.h_), variables_(ad.variables_), index_(ad.index_), der1_(ad.der1_), der2_(ad.der2_), crossDer2_(ad.crossDer2_),
      computeD1_(ad.computeD1_), computeD2_(ad.computeD2_), computeCrossD2_(ad.computeCrossD2_) {}

    AbstractNumericalDerivative& operator=(const AbstractNumericalDerivative& ad)
    {
      FunctionWrapper::operator=(ad);
      function1_ = ad.function1_;
      function2_ = ad.function2_;
      h_ = ad.h_;
      variables_ = ad.variables_;
      index_ = ad.index_;
      der1_ = ad.der1_;
      der2_ = ad.der2_;
      crossDer2_ = ad.crossDer2_;
      computeD1_ = ad.computeD1_;
      computeD2_ = ad.computeD2_;
      computeCrossD2_ = ad.computeCrossD2_;
      return *this;
    }

    virtual ~AbstractNumericalDerivative() {}

    AbstractNumericalDerivative* clone() const = 0;

  public:
    /**
     * @brief Set the interval value used in numerical approximation.
     *
     * Default value is 0.0001.
     *
     * @param h Interval value.
     */
    void setInterval(double h) { h_ = h; }
    
    /**
     * @return The interval value used in numerical approximation.
     */
    double getInterval() const { return h_; }
    
    /**
     * @brief Set the list of parameters to derivate numerically.
     *
     * @param variables A list of all parameter names.
     */
    void setParametersToDerivate(const std::vector<std::string>& variables)
    {
      variables_ = variables;
      index_.clear();
      for(size_t i = 0; i < variables_.size(); i++)
        index_[variables_[i]] = i;
      der1_.resize(variables_.size());
      der2_.resize(variables_.size());
      crossDer2_.resize(variables_.size(), variables_.size());
    }
    
    /**
     * @name The DerivableFirstOrder interface
     *
     * @{
     */
    void enableFirstOrderDerivatives(bool yn) { computeD1_ = yn; }
    bool enableFirstOrderDerivatives() const { return computeD1_; }
    
    double getFirstOrderDerivative(const std::string& variable) const
      throw (Exception)
    {
      std::map<std::string, size_t>::iterator it = index_.find(variable);
      if (computeD1_ && it != index_.end())
        return der1_[it->second];

      if (function1_)
        return function1_->getFirstOrderDerivative(variable);

      throw Exception("First order derivative not computed for variable " + variable + "."); 
    }
    /** @} */
    
    /**
     * @name The DerivableSecondOrder interface
     *
     * @{
     */

    void enableSecondOrderDerivatives(bool yn) { computeD2_ = yn; }
    bool enableSecondOrderDerivatives() const { return computeD2_; }

    double getSecondOrderDerivative(const std::string& variable) const
      throw (Exception)
    {
      std::map<std::string, size_t>::iterator it = index_.find(variable);
      if(computeD2_ && it != index_.end())
        return der2_[it->second];

      if (function2_)
        return function2_->getSecondOrderDerivative(variable);

      throw Exception("Second order derivative not computed for variable " + variable + "."); 
    }

    double getSecondOrderDerivative(const std::string& variable1, const std::string& variable2) const
      throw (Exception)
    {
      std::map<std::string, size_t>::iterator it1 = index_.find(variable1);
      std::map<std::string, size_t>::iterator it2 = index_.find(variable2);
      if(computeCrossD2_ && it1 != index_.end() && it2 != index_.end()) return crossDer2_(it1->second, it2->second);

      if (function2_)
        return function2_->getSecondOrderDerivative(variable1, variable2);

      throw Exception("Cross second order derivative not computed for variables " + variable1 + " and " + variable2 + "."); 
    }
  
    /** @} */
     
    /**
     * @name The Parametrizable interface.
     *
     * @{
     */
    double f(const ParameterList& parameters) throw (Exception)
    {
      setParameters(parameters);
      return getValue();
    }
    void setParameters(const ParameterList& parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
      function_->setParameters(parameters);
      updateDerivatives(parameters);
    }
    void setAllParametersValues(const ParameterList& parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
      function_->setAllParametersValues(parameters);
      updateDerivatives(parameters);
    }
    
    void setParameterValue(const std::string& name, double value)
      throw (ParameterNotFoundException, ConstraintException)
    {
      function_->setParameterValue(name, value);
      updateDerivatives(function_->getParameters().subList(name));
    }
    
    void setParametersValues(const ParameterList& parameters)
      throw (ParameterNotFoundException, ConstraintException)
    {
      function_->setParametersValues(parameters);
      updateDerivatives(parameters);
    }
    
    bool matchParametersValues(const ParameterList& parameters)
      throw (ConstraintException)
    {
      bool test = function_->matchParametersValues(parameters);
      updateDerivatives(parameters);
      return test;
    }
    /** @} */

    void enableSecondOrderCrossDerivatives(bool yn) { computeCrossD2_ = yn; }
    bool enableSecondOrderCrossDerivatives() const { return computeCrossD2_; }

  protected:
    /**
     * @brief Compute derivatives.
     *
     * @param parameters The point where to compute derivatives. It is NOT passed as references,
     * as the inner parameters of the function will be changed when computing the numerical derivatives.
     */
    virtual void updateDerivatives(const ParameterList parameters) = 0;
    
};

} //end of namespace bpp.

#endif //_ABSTRACTNUMERICALDERIVATIVE_H_

