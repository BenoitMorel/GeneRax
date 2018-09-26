//
// File: FunctionOperator.h
// Created by: Laurent Guéguen
// Created on: lundi 5 décembre 2016, à 23h 05
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

#ifndef _FUNCTION_OPERATOR_H_
#define _FUNCTION_OPERATOR_H_

#include <memory>
#include "../Functions.h"

namespace bpp
{

/**
 * @brief Implements a double operator (ie leaf in the computation
 * tree) where value comes from a function.
 *
 */

  template<class F>
  class FunctionOperator:
    public Operator
  {
  private:
    F& func_;

    std::string name_;

    double getValue_(std::true_type) const 
    {
      return func_.getValue();
    }

    double getValue_(std::false_type) const
    {
      return 0;
    }

    double getFirstOrderDerivative_(const std::string& variable, std::true_type) const
    {
      return func_.getFirstOrderDerivative(variable);
    }

    double getFirstOrderDerivative_(const std::string& variable, std::false_type) const
    {
      return 0;
    }

    double getSecondOrderDerivative_(const std::string& variable, std::true_type) const
    {
      return func_.getSecondOrderDerivative(variable);
    }

    double getSecondOrderDerivative_(const std::string& variable, std::false_type) const
    {
      return 0;
    }


  public:

    FunctionOperator(F& func, std::string name) :
      func_(func), name_(name)
    {
    }
                                
    FunctionOperator* clone() const 
    {
      return new FunctionOperator(*this);
    }

    double getValue() const
    {
      return getValue_(std::integral_constant<bool, std::is_base_of<Function, F>::value>{});
    }

    double getFirstOrderDerivative(const std::string& variable) const
    {
      return getFirstOrderDerivative_(variable, std::integral_constant<bool, std::is_base_of<DerivableFirstOrder, F>::value>{});
    }

    double getSecondOrderDerivative(const std::string& variable) const
    {
      return getSecondOrderDerivative_(variable, std::integral_constant<bool, std::is_base_of<DerivableSecondOrder, F>::value>{});
    }

    std::string output() const
    {
      return name_;
    }
    
  };


} //end of namespace bpp.

#endif  //_FUNCTION_OPERATOR_H_

