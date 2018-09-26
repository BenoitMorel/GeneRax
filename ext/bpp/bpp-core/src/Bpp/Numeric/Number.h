//
// File: Number.h
// Created by: Julien Dutheil
// Created on: Thu Nov 13 16:29:03 2003
//


/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide utilitary
classes. This file belongs to the Bio++ Project.

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

#ifndef _NUMBER_H_
#define _NUMBER_H_

#include "../Clonable.h"
#include "../Text/TextTools.h"

#include <string>

namespace bpp
{
/**
 * @brief The Number interface.
 *
 * This template class may be used to deal with number in an object way.
 */
class BppNumberI: public Clonable
{
	public:
		
		BppNumberI() {}
			
		virtual ~BppNumberI() {}

	public:
	
		virtual BppNumberI* clone() const = 0;
	
	public:
		
    virtual std::string toString() const = 0;

};


class BppNotANumber: public virtual BppNumberI
{
	public:
		
		BppNotANumber() {}
			
		virtual ~BppNotANumber() {}

	public:
	
		virtual BppNotANumber* clone() const { return new BppNotANumber(); }
	
	public:
		
    virtual std::string toString() const { return "NaN"; }

};


/**
 * @brief The Number object template class.
 *
 * This template class may be used to deal with number in an object way.
 */
template<class T> class Number: public virtual BppNumberI
{
	protected:
		/** @brief The value of this parameter. */
		T value_;
	
	public:
		
		/**
		 * @brief Build a new Number object with a specific value.
		 *
		 * @param value The value that the Number must have.
		 */
		Number(const T& value = 0): value_(value) {}
			
		virtual ~Number() {}

    Number<T> & operator=(const T & t)
    { 
      value_ = t;
      return *this;
    }
	
	public:
	
		/**
		 * @name The Clonable interface.
		 *
		 * @{
		 */
		Number<T>* clone() const { return new Number<T>(value_); }
		/** @} */
	
	public:
		
		/**
		 * @brief Get the value of this number.
		 *
		 * @return The value of this number.
		 */
		T getValue() const { return value_; }

    std::string toString() const { return TextTools::toString(value_); }
};

/**
 * @brief An object wrapper for double values.
 */
class BppDouble: public virtual Number<double>
{
 	public:
		
		/**
		 * @brief Build a new BppDouble number object with a specific value.
		 *
		 * @param value The value that the Number must have.
		 */
		BppDouble(double value = 0): Number<double>(value) {}
			
		virtual ~BppDouble() {}

	public:
	
		/**
		 * @name The Clonable interface.
		 *
		 * @{
		 */
		BppDouble* clone() const { return new BppDouble(*this); }
		/** @} */
	
};

/**
 * @brief An object wrapper for integer values.
 */
class BppInteger: public virtual Number<int>
{
 	public:
		
		/**
		 * @brief Build a new BppInteger number object with a specific value.
		 *
		 * @param value The value that the Number must have.
		 */
		BppInteger(int value = 0): Number<int>(value) {}
			
		virtual ~BppInteger() {}

	public:
	
		/**
		 * @name The Clonable interface.
		 *
		 * @{
		 */
		BppInteger* clone() const { return new BppInteger(*this); }
		/** @} */
	
};

/**
 * @brief An object wrapper for unsigned integer values.
 */
class BppUnsignedInteger: public virtual Number<unsigned int>
{
 	public:
		
		/**
		 * @brief Build a new BppUnsignedInteger number object with a specific value.
		 *
		 * @param value The value that the Number must have.
		 */
		BppUnsignedInteger(unsigned int value = 0): Number<unsigned int>(value) {}
			
		virtual ~BppUnsignedInteger() {}

	public:
	
		/**
		 * @name The Clonable interface.
		 *
		 * @{
		 */
		BppUnsignedInteger* clone() const { return new BppUnsignedInteger(*this); }
		/** @} */
	
};

} //end of namespace bpp.

#endif	//_NUMBER_H_

