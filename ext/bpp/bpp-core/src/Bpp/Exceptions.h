//
// File Exceptions.h
// Created by: Guillaume Deuchst
//              Julien Dutheil
//              Sylvain Gaillard
// Last modification : Thu Jul 22 2004
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

#ifndef _EXCEPTIONS_H_
#define _EXCEPTIONS_H_

#include <stdexcept>
#include <vector>
#include <string>
#include "Text/TextTools.h"

namespace bpp
{

/**
 * @brief Exception base class.
 *
 * Overload exception constructor (to control the exceptions mechanism).</p>
 */
class Exception:
  public std::exception
{
  protected:
    std::string message_;
  
  public:
    /**
     * @brief Build a new Exception.
     *
     * @param text A message to be passed to the exception hierarchy.
     */
    Exception(const char* text): message_(std::string(text)) {}

    /**
     * @brief Build a new Exception.
     *
     * @param text A message to be passed to the exception hierarchy.
     */
    Exception(const std::string& text): message_(text) {}
  
    virtual ~Exception() throw() {}

  public:
    
    /**
     * @brief Method to get the message of the exception (STL method redefinition).
     *
     * @return The message passed to the exception hierarchy.
     */
    const char* what() const throw() { return message_.c_str(); }
};


/**
 * @brief The base class exception for IO error.
 */
class IOException:
  public Exception
{
  public: // Class constructors and destructor:
    
    /**
     * @brief Build a new IOException.
     *
     * @param text A message to be passed to the exception hierarchy.
     */
    IOException(const char* text): Exception(text) {}

    /**
     * @brief Build a new IOException.
     *
     * @param text A message to be passed to the exception hierarchy.
     */
    IOException(const std::string& text): Exception(text) {} 
  
    virtual ~IOException() throw() {}

};


/**
 * @brief The base class exception for NULL pointer error.
 *
 * This exception may be thrown when an unexpected NULL pointer is found.
 */
class NullPointerException:
  public Exception
{
  public:

    /**
     * @brief Build a new NullPointerException.
     *
     * @param text A message to be passed to the exception hierarchy.
     */
    NullPointerException(const char* text): Exception(text) {}

    /**
     * @brief Build a new NullPointerException.
     *
     * @param text A message to be passed to the exception hierarchy.
     */
    NullPointerException(const std::string& text): Exception(text) {}

    virtual ~NullPointerException() throw() {}
};


/**
 * @brief The base class exception for zero division error.
 */
class ZeroDivisionException:
  public Exception
{
  public:

    /**
     * @brief Build a new ZeroDivisionException.
     *
     * @param text A message to be passed to the exception hierarchy.
     */
    ZeroDivisionException(const char* text): Exception(text) {}

    /**
     * @brief Build a new ZeroDivisionException.
     *
     * @param text A message to be passed to the exception hierarchy.
     */
    ZeroDivisionException(const std::string& text): Exception(text) {}

    virtual ~ZeroDivisionException() throw() {}
};

  
/**
 * @brief Number exception: integers.
 */
class BadIntegerException:
  public Exception
{
  protected:
    int badInt_;
  
  public:
    
    /**
     * @brief Build a new BadIntegerException.
     *
     * @param text   A message to be passed to the exception hierarchy.
     * @param badInt The faulty integer.
     */
    BadIntegerException(const char* text, int badInt);

    /**
     * @brief Build a new BadIntegerException.
     *
     * @param text   A message to be passed to the exception hierarchy.
     * @param badInt The faulty integer.
     */
    BadIntegerException(const std::string& text, int badInt);
  
    virtual ~BadIntegerException() throw() {}

  public:
    
    /**
     * @brief Get the integer that threw this exception.
     *
     * @return The faulty integer.
     */
    int getBadInteger() const { return badInt_; }
  
};


/**
 * @brief Number exception: doubles.
 */
class BadNumberException:
  public Exception
{
  protected:
    double badNumber_;
  
  public:
    
    /**
     * @brief Build a new BadNumberException.
     *
     * @param text      A message to be passed to the exception hierarchy.
     * @param badNumber The faulty number.
     */
    BadNumberException(const char* text, double badNumber);

    /**
     * @brief Build a new BadNumberException.
     *
     * @param text      A message to be passed to the exception hierarchy.
     * @param badNumber The faulty number.
     */
    BadNumberException(const std::string& text, double badNumber);
  
    virtual ~BadNumberException() throw() {}

  public:
    
    /**
     * @brief Get the number that threw this exception.
     *
     * @return The faulty number.
     */
    double getBadNumber() const { return badNumber_; }
  
};


/**
 * @brief Number format exception.
 */
class NumberFormatException:
  public Exception
{
  protected:
    std::string badNumber_;
  
  public:
    
    /**
     * @brief Build a new NumberFormatException.
     *
     * @param text      A message to be passed to the exception hierarchy.
     * @param badNumber The faulty number.
     */
    NumberFormatException(const char* text, const std::string& badNumber);

    /**
     * @brief Build a new NumberFormatException.
     *
     * @param text      A message to be passed to the exception hierarchy.
     * @param badNumber The faulty number.
     */
    NumberFormatException(const std::string& text, const std::string& badNumber);
  
    virtual ~NumberFormatException() throw() {}

  public:
    
    /**
     * @brief Get the number that threw this exception.
     *
     * @return The faulty number.
     */
    std::string getBadNumber() const { return badNumber_; }
  
};


/**
 * @brief Index out of bounds exception class.
 */
class IndexOutOfBoundsException:
  public virtual Exception
{
  protected:
    size_t badIndex_, lowerBound_, upperBound_;
  
  public:

    /**
     * @brief Build a new IndexOutOfBoundsException.
     *
     * @param text   A message to be passed to the exception hierarchy.
     * @param badInt The faulty integer.
     * @param lowerBound Lower limit.
     * @param upperBound Upper limit.
     */
  IndexOutOfBoundsException(const std::string& text, size_t badInt, size_t lowerBound, size_t upperBound) :
    Exception("out of [" + TextTools::toString(lowerBound) + ", " + TextTools::toString(upperBound) +  "])" + text),
    badIndex_(badInt),
    lowerBound_(lowerBound),
    upperBound_(upperBound) {}
  
    virtual ~IndexOutOfBoundsException() throw() {}

  public:
    
    /**
     * @brief Get the bounds.
     *
     * @return The bounds.
     */
    std::vector<size_t> getBounds() const;

    size_t getBadIndex() const { return badIndex_; }
};


/**
 * @brief Wrong size exception class.
 */
class BadSizeException:
  public virtual Exception
{
  protected:
    size_t badSize_, correctSize_;
  
  public:
    
    /**
     * @brief Build a new BadSizeException.
     *
     * @param text   A message to be passed to the exception hierarchy.
     * @param badSize The faulty size.
     * @param correctSize The expected size.
     */
    BadSizeException(const std::string& text, size_t badSize, size_t correctSize);
  
    virtual ~BadSizeException() throw() {}

  public:
    
    size_t getBadSize() const { return badSize_; }
    size_t getCorrectSize() const { return correctSize_; }
};


/**
 * @brief Out of range exception class.
 */
class OutOfRangeException:
  public Exception
{
  protected:
    double lowerBound_, upperBound_;
  
  public:
    
    /**
     * @brief Build a new OutOfRangeException.
     *
     * @param text   A message to be passed to the exception hierarchy.
     * @param badValue The faulty value.
     * @param lowerBound Lower limit.
     * @param upperBound Upper limit.
     */
    OutOfRangeException(const std::string& text, double badValue, double lowerBound, double upperBound);
  
    virtual ~OutOfRangeException() throw() {}

  public:
    
    /**
     * @return The lower bound.
     */
    double getLowerBound() const { return lowerBound_; }

    /**
     * @return The upper bound.
     */
    double getUpperBound() const { return upperBound_; }
};


/**
 * @brief This expeption is sent when a given method is not implemented.
 */
class NotImplementedException:
  public Exception
{
  public:

    /**
     * @brief Build a new NotImplementedException.
     *
     * @param text A message to be passed to the exception hierarchy.
     */
    NotImplementedException(const std::string& text): Exception(text) {}

    virtual ~NotImplementedException() throw() {}
};



} //end of namespace bpp.

#endif // _EXCEPTIONS_H_

