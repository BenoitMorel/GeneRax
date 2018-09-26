//
// File: Constraints.h
// Created by: Julien Dutheil
// Created on: Thu Dec 25 19:35:17 2003
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

#ifndef _CONSTRAINTS_H_
#define _CONSTRAINTS_H_

// From the STL:
#include <string>
#include <iostream>
#include <typeinfo>

// From Utils:
#include "../Clonable.h"
#include "../Text/TextTools.h"
#include "../Exceptions.h"

#include "NumConstants.h"

namespace bpp
{
/**
 * @brief The constraint interface.
 *
 * It provides a method that tells if a given value is correct.
 */
class Constraint : public Clonable
{
public:
  Constraint() {}
  virtual ~Constraint() {}

  Constraint* clone() const = 0;

public:
  /**
   * @brief Tell if a given value is correct.
   *
   * @param value The value to test.
   * @return True is the value is correct.
   */
  virtual bool isCorrect(double value) const = 0;

  /**
   * @brief Tell if all the values in a given interval are correct.
   *
   * @param min, max  The bounds of the interval.
   * @return True is the value is correct.
   */
  virtual bool includes(double min, double max) const = 0;

  /**
   * @brief Give the nearest limit for a bad value.
   *
   * @param value The bad value.
   * @return The nearer limit.
   */
  virtual double getLimit(double value) const = 0;

  /**
   * @brief Give the nearest accepted limit for a bad value.
   *
   * The difference with getLimit() is when the Constraint is open at
   * the limit, in which case the retruned value is the limit +- 1e-12.
   *
   * @param value The bad value.
   * @return The nearer limit.
   */
  virtual double getAcceptedLimit(double value) const = 0;

  /**
   * @brief Give a short description on the type of constraint.
   *
   * @return A string which describes the constraint.
   */
  virtual std::string getDescription() const = 0;

  /**
   * @brief Intersect this Constraint with another one
   *
   * @param c the intersected Constraint
   * @return the intersection
   */
  virtual Constraint* operator&(const Constraint& c) const = 0;

  /**
   * @brief Tells if this constraints defines an empty set
   */
  
  virtual bool isEmpty() const = 0;
};

/**
 * @brief An interval, either bounded or not, which can also have infinite bounds.
 *
 * The upper and lower bound can be included or not (strict bound), finite or infinite (in that case, equal to a very large value).
 * Despite the mathematical non-sense, and infinite bound can be either excluded or included.
 */

class IntervalConstraint : public Constraint
{
protected:
  /**
   * @brief The boundaries of the interval
   *
   **/
  double lowerBound_, upperBound_;

  /**
   * @brief Boolean flags are true if the boundaries are included
   *
   **/
  bool inclLowerBound_, inclUpperBound_;
  /**
   *
   * @brief the accepted precision on the boundary (default: 1e-12)
   **/

  double precision_;
  

public:
  IntervalConstraint() :  lowerBound_(NumConstants::MINF()),
                upperBound_(NumConstants::PINF()),
                inclLowerBound_(true),
                inclUpperBound_(true),
                precision_(NumConstants::TINY()) {}

  IntervalConstraint(double lowerBound, double upperBound, bool inclLower, bool inclUpper, double precision = NumConstants::TINY()) :
    lowerBound_(lowerBound),
    upperBound_(upperBound),
    inclLowerBound_(inclLower),
    inclUpperBound_(inclUpper),
    precision_(precision) {}
 
  /**
   * @brief Create an interval with an infinite lower/upper bound.
   *
   * The infinite bound will not be included, following mathematical conventions.
   *
   * @param isPositive Tell if the infinite bound is positive or negative.
   * @param bound The finite bound.
   * @param incl Tell if the finite bound is included or not.
   * @param precision Parameter precision.
   */

  IntervalConstraint(bool isPositive, double bound, bool incl, double precision = NumConstants::TINY()) :
    lowerBound_(isPositive ? bound : NumConstants::MINF()),
    upperBound_(isPositive ? NumConstants::PINF() : bound),
    inclLowerBound_(isPositive ? incl : false),
    inclUpperBound_(isPositive ? false : incl),
    precision_(precision) {}

  /**
   * @brief Create an interval from a string description, using
   * readDescription method.
   *
   **/

  IntervalConstraint(std::string& desc) :
    lowerBound_(NumConstants::MINF()),
    upperBound_(NumConstants::PINF()),
    inclLowerBound_(true),
    inclUpperBound_(true),
    precision_(NumConstants::TINY())
  {
    readDescription(desc);
  }
  
  virtual ~IntervalConstraint() {}

  IntervalConstraint* clone() const { return new IntervalConstraint(*this); }

public:
  void setLowerBound(double lowerBound, bool strict) { lowerBound_ = lowerBound; inclLowerBound_ = !strict; }
  void setUpperBound(double upperBound, bool strict) { upperBound_ = upperBound; inclUpperBound_ = !strict; }

  double getLowerBound() const { return lowerBound_; }
  double getUpperBound() const { return upperBound_; }

  bool strictLowerBound() const { return !inclLowerBound_; }
  bool strictUpperBound() const { return !inclUpperBound_; }

  bool finiteLowerBound() const { return lowerBound_ > NumConstants::MINF(); }
  bool finiteUpperBound() const { return upperBound_ < NumConstants::PINF(); }

  bool includes(double min, double max) const
  {
    return (inclLowerBound_ ? min >= getLowerBound() : min > getLowerBound()) &&
           (inclUpperBound_ ? max <= getUpperBound() : max < getUpperBound());
  }

  virtual bool isCorrect(double value) const
  {
    return (inclLowerBound_ ? value >= getLowerBound() : value > getLowerBound()) &&
           (inclUpperBound_ ? value <= getUpperBound() : value < getUpperBound());
  }

  bool operator<(double value) const
  {
    return inclUpperBound_ ? upperBound_ < value : upperBound_ <= value;
  }

  bool operator>(double value) const
  {
    return inclLowerBound_ ? lowerBound_ > value : lowerBound_ >= value;
  }

  bool operator<=(double value) const
  {
    return upperBound_ <= value;
  }

  bool operator>=(double value) const
  {
    return lowerBound_ >= value;
  }

  double getLimit(double value) const
  {
    return isCorrect(value) ? value :
           (*this >= value ? lowerBound_ : upperBound_);
  }

  double getAcceptedLimit(double value) const
  {
    return isCorrect(value) ? value :
           (*this >= value ?
            strictLowerBound() ? lowerBound_ + precision_ : lowerBound_ :
            strictUpperBound() ? upperBound_ - precision_ : upperBound_);
  }

  double getPrecision() const
  {
    return precision_;
  }
  
  std::string getDescription() const
  {
    return (inclLowerBound_ ? "[ " : "]")
           + (finiteLowerBound() ? TextTools::toString(lowerBound_) : "-inf")
           + "; "
           + (finiteUpperBound() ? TextTools::toString(upperBound_) : "+inf")
           + (inclUpperBound_ ? "] " : "[");
  }

  /**
   *
   * @brief Sets the bounds of the interval from a string.
   *
   * @param desc the description in interval-like syntax, with signs
   * "[", ";", "]" as well as floats and "-inf" and "inf".
   *
   **/
  
  void readDescription(std::string& desc)
  {
    size_t pdp=desc.find(";");
    size_t dc=desc.find_first_of("[]",1);

    if (dc==std::string::npos || pdp==std::string::npos ||
        (desc[0]!=']' && desc[0]!='[') || (pdp >=dc))
      throw Exception("Constraints::readDescription. Wrong description:" + desc);

    std::string deb=desc.substr(1,pdp-1);
    std::string fin=desc.substr(pdp+1, dc-pdp-1);

    inclLowerBound_=(desc[0]=='[');
    inclUpperBound_=(desc[dc]==']');

    lowerBound_ =  (deb=="-inf")?NumConstants::MINF():TextTools::toDouble(deb);
    upperBound_ =  ((fin=="+inf") || (fin=="inf"))?NumConstants::PINF():TextTools::toDouble(fin);
    
  }
  
  /**
   * @brief Intersect this IntervalConstraint with another one
   *
   * @param c the intersected IntervalConstraint
   * @return the intersection, or NULL if c is not an IntervalConstraint. The
   * resulting precision is the maximum of both precisions.
   */
  Constraint* operator&(const Constraint& c) const
  {
    double lowerBound, upperBound;
    bool inclLowerBound, inclUpperBound;

    const IntervalConstraint* pi = dynamic_cast<const IntervalConstraint*>(&c);

    if (pi)
    {
      if (lowerBound_ <= pi->lowerBound_)
      {
        lowerBound = pi->lowerBound_;
        inclLowerBound = pi->inclLowerBound_;
      }
      else
      {
        lowerBound = lowerBound_;
        inclLowerBound = inclLowerBound_;
      }

      if (upperBound_ >= pi->upperBound_)
      {
        upperBound = pi->upperBound_;
        inclUpperBound = pi->inclUpperBound_;
      }
      else
      {
        upperBound = upperBound_;
        inclUpperBound = inclUpperBound_;
      }
      return new IntervalConstraint(lowerBound, upperBound, inclLowerBound, inclUpperBound, (precision_>pi->getPrecision())?precision_:pi->getPrecision());
    }
    else
      return 0;
  }

  /**
   * @brief Intersect this IntervalConstraint with another constraint
   *
   * @param c the intersected constraint
   * @return this IntervalConstraint modified, or not modified if c is not an
   * IntervalConstraint. The precision is set to the maximum of bith precisions.
   */
  IntervalConstraint& operator&=(const Constraint& c)
  {
    try {
      const IntervalConstraint& pi = dynamic_cast<const IntervalConstraint&>(c);

      if (lowerBound_ <= pi.lowerBound_)
      {
        lowerBound_ = pi.lowerBound_;
        inclLowerBound_ = pi.inclLowerBound_;
      }

      if (upperBound_ >= pi.upperBound_)
      {
        upperBound_ = pi.upperBound_;
        inclUpperBound_ = pi.inclUpperBound_;
      }
      if (pi.getPrecision() > precision_)
        precision_ = pi.getPrecision();
    } catch(std::bad_cast&) {}

    return *this;
  }

  /**
   * @brief Tells if this interval equals another one
   *
   * @param i the compared IntervalConstraint
   */
  bool operator==(const IntervalConstraint& i) const
  {
    return lowerBound_ == i.lowerBound_
      && inclLowerBound_ == i.inclLowerBound_
      && upperBound_ == i.upperBound_
      && inclUpperBound_ == i.inclUpperBound_;
  }

  /**
   * @brief Tells if this interval is different from another one
   *
   * @param i the compared IntervalConstraint
   */
  bool operator!=(const IntervalConstraint& i) const
  {
    return lowerBound_ != i.lowerBound_
      || inclLowerBound_ != i.inclLowerBound_
      || upperBound_ != i.upperBound_
      || inclUpperBound_ != i.inclUpperBound_;
  }

  /**
   * @brief Tells if this interval is included or equal in another one
   *
   * @param i the compared IntervalConstraint
   */
  bool operator<=(const IntervalConstraint& i) const
  {
    return lowerBound_ >= i.lowerBound_
           && upperBound_ <= i.upperBound_;
  }

  /**
   * @brief Tells if this interval is empty 
   */
  
  bool isEmpty() const
  {
    return ((lowerBound_ > upperBound_) ||
            ((lowerBound_ > upperBound_) &&
             inclUpperBound_ && inclLowerBound_));
  }

};


} // end of namespace bpp.

#endif  // _CONSTRAINTS_H_

