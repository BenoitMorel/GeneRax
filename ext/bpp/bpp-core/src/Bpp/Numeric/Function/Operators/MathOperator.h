//
// File: MathOperator.h
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

#ifndef _MATH_OPERATOR_H_
#define _MATH_OPERATOR_H_

#include <memory>
#include <cmath>

namespace bpp
{

/**
 * @brief Implements a unary operator that applies a math (described
 * in cmath) operator. 
 *
 */

  class MathOperator:
    public Operator
  {
  private:
    /*
     * @brief pointer to function.
     *  If null, identity function is used.
     */
    
    double (*func_)(double);

    std::string name_;

    std::shared_ptr<Operator> son_;
    
  public:

    MathOperator(double (*func)(double), std::string name, std::shared_ptr<Operator> son) :
      func_(func), name_(name), son_(son)
    {
    }

    MathOperator* clone() const 
    {
      return new MathOperator(*this);
    }
    

    double getValue() const
    {
      if (func_)
        return (*func_)(son_->getValue());
      else
        return son_->getValue();
    }

    /**
     * @brief 1st order derivative
     *
     */
    
    double getFirstOrderDerivative(const std::string& variable) const
    {
      double v = son_->getValue();
      double d = son_->getFirstOrderDerivative(variable);

      if (name_=="exp")
        return d*exp(v);
      else if (name_=="log")
        return d/v;
      else
        throw Exception("MathOperator::getFirstOrderDerivative : unknown function " + name_);
      
    }
    
    double getSecondOrderDerivative(const std::string& variable) const
    {
      double v = son_->getValue();
      double d = son_->getFirstOrderDerivative(variable);
      double d2 = son_->getSecondOrderDerivative(variable);

      if (name_=="exp")
        return (d2+d*d)*exp(v);
      else if (name_=="log")
        return (d2*v-d*d)/(v*v);
      else
        throw Exception("MathOperator::getFirstOrderDerivative : unknown function " + name_);
      return 0;
    }


    std::string output() const
    {
      return name_ + "(" + son_->output() + ")";
    }
    
    
  };
  

} //end of namespace bpp.

#endif  //_MATH_OPERATOR_H_

