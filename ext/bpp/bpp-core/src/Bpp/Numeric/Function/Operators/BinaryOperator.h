//
// File: BinaryOperator.h
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

#ifndef _BINARY_OPERATOR_H_
#define _BINARY_OPERATOR_H_

#include <memory>

namespace bpp
{

/**
 * @brief Binary arithmetic operator for numerical computation.
 *
 */

  class BinaryOperator:
    public Operator
  {
  private:
    char symb_;

    std::shared_ptr<Operator> left_, right_;
    
  public:

    BinaryOperator(char symb, std::shared_ptr<Operator> left, std::shared_ptr<Operator> right) :
      symb_(symb), left_(left), right_(right)
    {
    }

    BinaryOperator* clone() const 
    {
      return new BinaryOperator(*this);
    }

    char getSymbol() const
    {
      return symb_;
    }
    
    double getValue() const
    {
      switch(symb_)
      {
      case '+':
        return left_->getValue() + right_->getValue();
      case '-':
        return left_->getValue() - right_->getValue();
      case '/':
        if (right_->getValue()==0)
          return 0;
        
        return left_->getValue() / right_->getValue();
      case '*':
        return left_->getValue() * right_->getValue();
      default:
        return 0;
      }
    }

    double getFirstOrderDerivative(const std::string& variable) const
    {
      double dl = left_->getFirstOrderDerivative(variable);
      double dr = right_->getFirstOrderDerivative(variable);
      double l = left_->getValue();
      double r = right_->getValue();

      switch(symb_)
      {
      case '+':
        return dl + dr;
      case '-':
        return dl - dl;
      case '/':
        if (r==0)
          return 0;
        
        return (dl * r - dr * l ) / (r * r);
      case '*':
        return dl * r + dr * l;
      default:
        return 0;
      }
      return 0;
    }
    
    double getSecondOrderDerivative(const std::string& variable) const
    {
      double d2l = left_->getSecondOrderDerivative(variable);
      double d2r = right_->getSecondOrderDerivative(variable);
      double l = left_->getValue();
      double r = right_->getValue();
      double dl = left_->getFirstOrderDerivative(variable);
      double dr = right_->getFirstOrderDerivative(variable);

      double r2 = r*r;
      double r3 = r*r2;
      
      switch(symb_)
      {
      case '+':
        return d2l + d2r;
      case '-':
        return d2l - d2r;
        
      case '/':
        if (r==0)
          return 0;
        
        return (d2l * r - d2r * l ) / r2 - ( 2 * dl * ( dl * r - dr * l) ) / r3;
        
      case '*':
        return d2l * r + d2r * l + 2 * dr * dl;
      default:
        return 0;
      }
    }

    std::string output() const
    {
      return "(" + left_->output() + " " + symb_ + " " + right_->output() + ")";
    }
    

  };
  

} //end of namespace bpp.

#endif  //_BINARY_OPERATOR_H_

