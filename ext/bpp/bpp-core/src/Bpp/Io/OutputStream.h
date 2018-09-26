//
// File: OutputStream.h
// Created by: Julien Dutheil
// Created on: Mon Jan 25 17:41 2010
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

#ifndef _OUTPUTSTREAM_H_
#define _OUTPUTSTREAM_H_

#include "../Clonable.h"

//From the STL:
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <memory>
#include <iomanip>

namespace bpp
{

/**
 * @brief OutputStream interface.
 *
 * This interface define several << operators for outputing messages.
 * Wrapper classes for the STL streams are available as special
 * implementations, but this interface can be used to redirect the output to
 * files or to graphical components for instance.
 */
class OutputStream :
  public virtual Clonable
{
public:
  virtual OutputStream& operator<<(const std::string& message) = 0;
  virtual OutputStream& operator<<(const char* message) = 0;
  virtual OutputStream& operator<<(const char& message) = 0;
  virtual OutputStream& operator<<(const int& message) = 0;
  virtual OutputStream& operator<<(const unsigned int& message) = 0;
  virtual OutputStream& operator<<(const long int& message) = 0;
  virtual OutputStream& operator<<(const unsigned long int& message) = 0;
  virtual OutputStream& operator<<(const double& message) = 0;
  virtual OutputStream& operator<<(const long double& message) = 0;
  virtual OutputStream& operator<<(const bool& message) = 0;
  virtual OutputStream& endLine() = 0;
  virtual OutputStream& flush() = 0;
  virtual OutputStream& setPrecision(int digit) = 0;
  virtual int getPrecision() const = 0;
  virtual OutputStream& enableScientificNotation(bool yn) = 0;
  virtual bool isScientificNotationEnabled() const = 0;

  /**
   * @name The Clonable interface.
   *
   * @{
   */
  OutputStream* clone() const = 0;
  /** @} */

};





/**
 * @brief Helper implementation of the OutputStream interface.
 */
class AbstractOutputStream :
  public virtual OutputStream
{
private:
  int precision_;
  bool scienceNotation_;

public:
  AbstractOutputStream() : precision_(6), scienceNotation_(false) {}

public:
  OutputStream& setPrecision(int digit)
  {
    precision_ = digit;
    return *this;
  }
  int getPrecision() const { return precision_; }

  virtual OutputStream& enableScientificNotation(bool yn) { scienceNotation_ = yn; return *this; }
  virtual bool isScientificNotationEnabled() const { return scienceNotation_; }
};





/**
 * @brief Null output stream (swich off output).
 */
class NullOutputStream :
  public AbstractOutputStream
{
public:
  NullOutputStream& operator<<(const std::string& message) { return *this; }
  NullOutputStream& operator<<(const char* message) { return *this; }
  NullOutputStream& operator<<(const char& message) { return *this; }
  NullOutputStream& operator<<(const int& message) { return *this; }
  NullOutputStream& operator<<(const unsigned int& message) { return *this; }
  NullOutputStream& operator<<(const long int& message) { return *this; }
  NullOutputStream& operator<<(const unsigned long int& message) { return *this; }
  NullOutputStream& operator<<(const double& message) { return *this; }
  NullOutputStream& operator<<(const long double& message) { return *this; }
  NullOutputStream& operator<<(const bool& message) { return *this; }
  NullOutputStream& endLine() { return *this; }
  NullOutputStream& flush() { return *this; }

  NullOutputStream* clone() const { return new NullOutputStream(*this); }

};






/**
 * @brief STL output stream.
 *
 * This class wraps the std::ostream class.
 * It takes as input a pointer toward an existing stream that will then be
 * owned by the wrapper, as a smart pointer.
 * Any copy of this class will then result in an inactivation of the original
 * version (in other word, you can't have two wrappers for the same stream).
 */
class StlOutputStream :
  public AbstractOutputStream
{
private:
  mutable std::unique_ptr<std::ostream> stream_;

public:
  StlOutputStream(std::ostream* stream): stream_(stream) {}
  StlOutputStream(const StlOutputStream& stlos) : stream_() { stream_ = std::move(stlos.stream_); }
  StlOutputStream& operator=(const StlOutputStream& stlos)
  {
    stream_ = std::move(stlos.stream_);
    return *this;
  }

public:
  StlOutputStream& operator<<(const std::string& message) { if (stream_.get()) *stream_ << message; return *this; }
  StlOutputStream& operator<<(const char* message) { if (stream_.get()) *stream_ << message; return *this; }
  StlOutputStream& operator<<(const char& message) { if (stream_.get()) *stream_ << message; return *this; }
  StlOutputStream& operator<<(const int& message) { if (stream_.get()) *stream_ << message; return *this; }
  StlOutputStream& operator<<(const unsigned int& message) { if (stream_.get()) *stream_ << message; return *this; }
  StlOutputStream& operator<<(const long int& message) { if (stream_.get()) *stream_ << message; return *this; }
  StlOutputStream& operator<<(const unsigned long int& message) { if (stream_.get()) *stream_ << message; return *this; }
  StlOutputStream& operator<<(const double& message)
  {
    if (stream_.get())
      *stream_ << std::setprecision(getPrecision()) << (isScientificNotationEnabled() ? std::scientific : std::fixed) << message;
    return *this;
  }
  StlOutputStream& operator<<(const long double& message)
  {
    if (stream_.get())
      *stream_ << std::setprecision(getPrecision()) << (isScientificNotationEnabled() ? std::scientific : std::fixed) << message;
    return *this;
  }
  StlOutputStream& operator<<(const bool& message) { if (stream_.get()) *stream_ << message; return *this; }
  StlOutputStream& endLine() { if (stream_.get()) *stream_ << std::endl; return *this; }
  StlOutputStream& flush() { if (stream_.get()) stream_->flush(); return *this; }

  StlOutputStream* clone() const { return new StlOutputStream(*this); }

};





/**
 * @brief STL wrapper for output stream.
 *
 * This class wraps the std::ostream class, by forwarding to the STL class.
 * It does not own the STL stream and won't delete it.
 */
class StlOutputStreamWrapper :
  public AbstractOutputStream
{
protected:
  std::ostream* stream_;

public:
  StlOutputStreamWrapper(std::ostream* stream): stream_(stream) {}
  StlOutputStreamWrapper(const StlOutputStreamWrapper& stlos) : stream_(stlos.stream_) {}
  StlOutputStreamWrapper& operator=(const StlOutputStreamWrapper& stlos) { stream_ = stlos.stream_; return *this; }

public:
  StlOutputStreamWrapper& operator<<(const std::string& message) { if (stream_) *stream_ << message; return *this; }
  StlOutputStreamWrapper& operator<<(const char* message) { if (stream_) *stream_ << message; return *this; }
  StlOutputStreamWrapper& operator<<(const char& message) { if (stream_) *stream_ << message; return *this; }
  StlOutputStreamWrapper& operator<<(const int& message) { if (stream_) *stream_ << message; return *this; }
  StlOutputStreamWrapper& operator<<(const unsigned int& message) { if (stream_) *stream_ << message; return *this; }

  StlOutputStreamWrapper& operator<<(const long int& message) { if (stream_) *stream_ << message; return *this; }
  StlOutputStreamWrapper& operator<<(const unsigned long int& message) { if (stream_) *stream_ << message; return *this; }
  StlOutputStreamWrapper& operator<<(const double& message)
  {
    if (stream_)
      *stream_ << std::setprecision(getPrecision()) << (isScientificNotationEnabled() ? std::scientific : std::fixed) << message;
    return *this;
  }
  StlOutputStreamWrapper& operator<<(const long double& message)
  {
    if (stream_)
      *stream_ << std::setprecision(getPrecision()) << (isScientificNotationEnabled() ? std::scientific : std::fixed) << message;
    return *this;
  } 
  StlOutputStreamWrapper& operator<<(const bool& message) { if (stream_) *stream_ << message; return *this; }
  StlOutputStreamWrapper& endLine() { if (stream_) *stream_ << std::endl; return *this; }
  StlOutputStreamWrapper& flush() { if (stream_) stream_->flush(); return *this; }

  StlOutputStreamWrapper* clone() const { return new StlOutputStreamWrapper(*this); }

};

/**
 * @brief Standard output stream.
 *
 * This class wraps the std::cout stream.
 */
class StdOut :
  public AbstractOutputStream
{
public:
  OutputStream& operator<<(const std::string& message) { std::cout << message; return *this; }
  OutputStream& operator<<(const char* message) { std::cout << message; return *this; }
  OutputStream& operator<<(const char& message) { std::cout << message; return *this; }
  OutputStream& operator<<(const int& message) { std::cout << message; return *this; }
  OutputStream& operator<<(const unsigned int& message) { std::cout << message; return *this; }
  OutputStream& operator<<(const long int& message) { std::cout << message; return *this; }
  OutputStream& operator<<(const unsigned long int& message) { std::cout << message; return *this; }
  OutputStream& operator<<(const double& message)
  {
    std::cout << std::setprecision(getPrecision()) << (isScientificNotationEnabled() ? std::scientific : std::fixed) << message;
    return *this;
  }
  OutputStream& operator<<(const long double& message)
  {
    std::cout << std::setprecision(getPrecision()) << (isScientificNotationEnabled() ? std::scientific : std::fixed) << message;
    return *this;
  }
  OutputStream& operator<<(const bool& message) { std::cout << message; return *this; }
  OutputStream& endLine() { std::cout << std::endl; return *this; }
  OutputStream& flush() { std::cout.flush(); return *this; }

  StdOut* clone() const { return new StdOut(*this); }

};





/**
 * @brief Standard error stream.
 *
 * This class wraps the std::cerr stream.
 */
class StdErr :
  public AbstractOutputStream
{
public:
  OutputStream& operator<<(const std::string& message) { std::cerr << message; return *this; }
  OutputStream& operator<<(const char* message) { std::cerr << message; return *this; }
  OutputStream& operator<<(const char& message) { std::cerr << std::setprecision(getPrecision()) << message; return *this; }
  OutputStream& operator<<(const int& message) { std::cerr << std::setprecision(getPrecision()) << message; return *this; }
  OutputStream& operator<<(const unsigned int& message) { std::cerr << message; return *this; }
  OutputStream& operator<<(const long int& message) { std::cerr << std::setprecision(getPrecision()) << message; return *this; }
  OutputStream& operator<<(const unsigned long int& message) { std::cerr << message; return *this; }
  OutputStream& operator<<(const double& message)
  {
    std::cerr << std::setprecision(getPrecision()) << (isScientificNotationEnabled() ? std::scientific : std::fixed) << message;
    return *this;
  }
  OutputStream& operator<<(const long double& message)
  {
    std::cerr << std::setprecision(getPrecision()) << (isScientificNotationEnabled() ? std::scientific : std::fixed) << message;
    return *this;
  }
  OutputStream& operator<<(const bool& message) { std::cerr << message; return *this; }
  OutputStream& endLine() { std::cerr << std::endl; return *this; }
  OutputStream& flush() { std::cerr.flush(); return *this; }

  StdErr* clone() const { return new StdErr(*this); }

};

/**
 * @brief String output stream.
 *
 * This class wraps the std::ostringstream stream.
 */
class StdStr :
  public StlOutputStreamWrapper
{
public:
  StdStr(): StlOutputStreamWrapper(new std::ostringstream()){}
  
  std::string str() const { return dynamic_cast<const std::ostringstream*>(stream_)->str();}

  ~StdStr() { delete stream_;}
};

} // end of namespace bpp;

#endif //_OUTPUTSTREAM_H_

