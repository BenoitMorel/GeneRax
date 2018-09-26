//
// File: TransformedParameter.h
// Created by: Julien Dutheil
// Created on: Fri Jan 30 09:42 2009
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

#ifndef _TRANSFORMEDPARAMETER_H_
#define _TRANSFORMEDPARAMETER_H_

#include "Parameter.h"
#include "NumConstants.h"

#include <cmath>

namespace bpp
{

/**
 * @brief The TransformedParameter abstract class.
 *
 * This class extends the Parameter class.
 * A transformed parameter does not have a constraint attached to it, and is supposed to range from -inf to +inf.
 * It uses a transformation in order to do this, typically using a bijection.
 * The exact function used to achieve the transformation depends on the implementation of the interface.
 */ 
class TransformedParameter:
  public Parameter
{
	public:
    TransformedParameter(const std::string& name, double value):
      Parameter(name, value) {}
		
		TransformedParameter* clone() const = 0;
	
	public:
    /**
     * @brief Set the value of the parameter using the orignal coordinate system.
     *
     * @param value Parameter value in original coordinates.
     * @throw ConstraintException if the value is not correct.
     */
    virtual void setOriginalValue(double value) throw (ConstraintException) = 0;

    /**
     * @return The current value of the parameter in orignal coordinates.
     */
    virtual double getOriginalValue() const = 0;

    /**
     * @return The first order derivative of the transformation at the original point.
     * @throw NotImplementedException if the transformation does not support derivation
     * or if the derivation was not implemented.
     */
    virtual double getFirstOrderDerivative() const throw (NotImplementedException) = 0;
    
    /**
     * @return The second order derivative of the transformation at the original point.
     * @throw NotImplementedException if the transformation does not support derivation
     * or if the derivation was not implemented.
     */
    virtual double getSecondOrderDerivative() const throw (NotImplementedException) = 0;
};

/**
 * @brief Parameter transformation from ] b, +inf [ or ] -inf, b [ to ]-inf, + inf [.
 *
 * The equation of the tranformation is 
 * @f[
 * x' = \begin{cases}
 *   \log(a\cdot(x-b)) & \text{if $x < b+1$},\\
 *   a(x-1-b)          & \text{if $a \geq b+1$}.
 * \end{cases}
 * @f]
 * for a transformation from ] b, +inf [ to ]-inf, + inf [.
 * The 'b' parameter is the lower bound and 'a' is a scaling factor set to 1 by default.
 * For a transformation from  ] -inf, b [, the transformation is then
 * @f[
 * x' = \begin{cases}
 *   -\log(-a\cdot(x-b)) & \text{if $x < b-1$},\\
 *   -a(x-1-b)           & \text{if $a \geq b-1$}.
 * \end{cases}
 * @f]
 */
class RTransformedParameter:
  public TransformedParameter
{
  private:
    double scale_;
    double bound_;
    bool positive_;

  public:
    /**
     * @brief Build a new RTransformedParameter, with given bound and scale.
     *
     * @param name the name of the parameter.
     * @param value the value of th eparameter, in orginal coordinates.
     * @param bound the inerval bound to use.
     * @param positive tell if the original interval is positive or negative.
     * @param scale the scaling factor.
     */
    RTransformedParameter(const std::string& name, double value, double bound = 0, bool positive = true, double scale = 1):
      TransformedParameter(name, 1.),
      scale_(scale),
      bound_(bound),
      positive_(positive)
    {
      setOriginalValue(value);
    }

    RTransformedParameter* clone() const { return new RTransformedParameter(*this); }

  public:
    void setOriginalValue(double value) throw (ConstraintException) 
    {
      if (positive_ ? value <= bound_ : value >= bound_) throw ConstraintException("RTransformedParameter::setValue", this, value);
      if (positive_  & (value < 1 + bound_)) setValue(log(scale_ * (value - bound_)));
      if (positive_  & (value >= 1 + bound_)) setValue(scale_ * (value - 1. - bound_));
      if (!positive_ & (value > -1 + bound_)) setValue(log(-scale_ * (value - bound_)));
      if (!positive_ & (value <= -1 + bound_)) setValue(-scale_ * (value - 1. - bound_));
    }

    double getOriginalValue() const
    {
      double x = getValue();
      if (positive_)
        if(x < 0) return exp(x) / scale_ + bound_;
        else      return x / scale_ + 1. + bound_;
      else 
        if(x < 0) return - exp(-x) / scale_ + bound_;
        else      return - x / scale_ - 1. + bound_;
    }

    double getFirstOrderDerivative() const throw (NotImplementedException)
    {
      double x = getValue();
      if (positive_)
        if(x < 0) return exp(x) / scale_;
        else      return 1. / scale_;
      else 
        if(x < 0) return exp(-x) / scale_;
        else      return - 1. / scale_;
    }

    double getSecondOrderDerivative() const throw (NotImplementedException)
    {
      double x = getValue();
      if (positive_)
        if(x < 0) return exp(x) / scale_;
        else      return 0;
      else 
        if(x < 0) return - exp(-x) / scale_;
        else      return 0;
    }

};

/**
 * @brief Parameter transformation from ] a, b [ to ]-inf, + inf [.
 *
 * The equation of the tranformation is 
 * @f[
 * x' = s\tan\left(\pi\frac{x-a}{b-a} - \frac{\pi}{2}\right)
 * @f]
 * The 'a' and 'b' parameters are the lower and upper bounds and 's' is a scaling factor set to 1 by default.
 * If the hyperbolic option is set to true (the default), then the following transformation is used instead:
 * @f[
 * x' = s\,\text{atanh}\left(2\frac{x-a}{b-a} - 1\right)
 * @f]
 *
 */
class IntervalTransformedParameter:
  public TransformedParameter
{
  private:
    double scale_;
    double lowerBound_;
    double upperBound_;
    bool hyper_;
    double tiny_;

  public:
    /**
     * @brief Build a new IntervalTransformedParameter, with given bounds and scale.
     *
     * @param name the name of the parameter.
     * @param value the value of th eparameter, in orginal coordinates.
     * @param lowerBound the inerval lower bound to use.
     * @param upperBound the inerval lower bound to use.
     * @param scale the scaling factor.
     * @param hyper tell if the hyberbolic function should be used (true by default).
     */
    IntervalTransformedParameter(const std::string& name, double value, double lowerBound = 0, double upperBound = 1, double scale = 1, bool hyper = true):
      TransformedParameter(name, hyper ?
          scale * atanh(2. * (value - lowerBound) / (upperBound - lowerBound) - 1.) :
          scale * tan(NumConstants::PI() * (value - lowerBound)/(upperBound - lowerBound) - NumConstants::PI() / 2.)),
      scale_(scale),
      lowerBound_(lowerBound),
      upperBound_(upperBound),
      hyper_(hyper),
      tiny_(NumConstants::TINY())
    {}

    IntervalTransformedParameter* clone() const { return new IntervalTransformedParameter(*this); }

  public:
    void setOriginalValue(double value) throw (ConstraintException) 
    {
      if (value <= lowerBound_ || value >= upperBound_) throw ConstraintException("IntervalTransformedParameter::setValue", this, value);
      setValue(hyper_ ?
          scale_ * atanh(2. * (value - lowerBound_) / (upperBound_ - lowerBound_) - 1.) :
          scale_ * std::tan(NumConstants::PI() * (value - lowerBound_)/(upperBound_ - lowerBound_) - NumConstants::PI() / 2.));
    }

    double getOriginalValue() const
    {
      double x = getValue();
      double x2 = hyper_ ?
        (tanh(x / scale_) + 1.) * (upperBound_ - lowerBound_) / 2. + lowerBound_ :
        (atan(x / scale_) + NumConstants::PI() / 2.) * (upperBound_ - lowerBound_) / NumConstants::PI() + lowerBound_;
      return x2;
    }


    double getFirstOrderDerivative() const throw (NotImplementedException)
    {
      double x = getValue();
      double x2 = hyper_ ?
        1. / (pow(cosh(x / scale_), 2)) * (upperBound_ - lowerBound_) / (2. * scale_) :
        (upperBound_ - lowerBound_) / (NumConstants::PI() * scale_ * (pow(x / scale_, 2) + 1.));
      return x2;
    }
    double getSecondOrderDerivative() const throw (NotImplementedException)
    {
      double x = getValue();
      double x2 = hyper_ ?
        - 1. / (pow(cosh(x / scale_), 2)) * tanh(x / scale_) * (upperBound_ - lowerBound_) / (scale_ * scale_) :
        -2. * x * (upperBound_ - lowerBound_) / (NumConstants::PI() * pow(scale_, 3) * pow((pow(x / scale_, 2) + 1.), 2));
      return x2;
    }
};

/**
 * @brief 'Placebo' parameter transformation from ] b, +inf [ or ] -inf, b [ to ]-inf, + inf [.
 *
 * The class create a Transformed parameter which is exactly the same as a standard parameter.
 * It only implements the setOriginalValue and getOriginalValue methods, and remove the constraint.
 */
class PlaceboTransformedParameter:
  public TransformedParameter
{
  public:
    PlaceboTransformedParameter(const std::string& name, double value):
      TransformedParameter(name, value)
    {}

    PlaceboTransformedParameter* clone() const { return new PlaceboTransformedParameter(*this); }

  public:
    void setOriginalValue(double value) throw (ConstraintException) 
    {
      setValue(value);
    }

    double getOriginalValue() const
    {
      return getValue();
    }
    
    double getFirstOrderDerivative() const throw (NotImplementedException) { return 1.; }
    
    double getSecondOrderDerivative() const throw (NotImplementedException) { return 0.; }
};


} //end of namespace bpp.

#endif	//_TRANSFORMEDPARAMETER_H_

