//
// File: Clonable.h
// Created by: Julien Dutheil
// Created on: Wed Nov 12 15:55:03 2003
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

#ifndef _CLONABLE_H_
#define _CLONABLE_H_

/**
 * @mainpage
 *
 * @par
 * This library contains Bio++ core interfaces, classes and functions.
 * It provides a general interface called bpp::Clonable, which allows dynamic copy of objects.
 * Wrapper classes for numbers, strings and vectors that implement this interface are provided, and called
 * respectively bpp::Number, bpp::BppString and bpp::BppVector.
 *
 * @par
 * The bpp::Exception class is the most general exception class used in Bio++, and provides only a single text information.
 * Several more specialized exception classes are also provided, like
 * - bpp::IOException for input/output errors, usually when reading from or writing to files,
 * - bpp::NullPointerException when a NULL pointer is unexpectedly provided,
 * - bpp::BadNumberException and bpp::BadIntegerException, dealing with number problems (bounds, format, etc),
 * - and more...
 *
 * @par
 * This library also includes several submodules:
 * - The Text module contains: functions to deal with character strings (bpp::StringTokenizer and the bpp::TextTools static class),
 * like pattern matching, concatenation, number conversion, etc.
 * - The Graphics module provides support for graphical output. This includes the general bpp::GraphicDevice interface, with currently three implementations for XFig (bpp::XFigGraphicDevice), SVG (bpp::SVGGraphicDevice) and PGF (bpp::PGFGraphicDevice) formats.
 * Classes to deal with colors and fonts were created (see bpp::RGBColor, bpp::ColorSet, bpp::ColorTools and bpp::Font).
 * - The App module contains classes and tools that provide powerful methods to parse command-line options and more.
 * - The Numeric module is quite large and provides interfaces, classes and functions for numerical calculus.
 *   The tools included are general enough to be useful beyond biology.
 *   The originality of this module, unless many existing libraries, is to be fully object-oriented.
 *   Available methods include:
 *   - Vector operations, based on the stl::vector class, see bpp::VectorTools,
 *   - Matrix data and operations, see the bpp::Matrix interface and bpp::MatrixTools class,
 *   - Functions are implemented using several general classes:
 *     - the bpp::Function interface and its derivatives,
 *     - the bpp::Parameter, bpp::Constraint and bpp::ParameterList classes, and the bpp::Parametrizable interface
 *       provide a very general way to deal with parameters
 *     - Function operations, see the bpp::Optimizer interface, and the bpp::AbstractNumericalDerivative class for
 *       numerical differenciation,
 *     - Discrete probability distributions are implemented via the bpp::DiscreteDistribution interface
 *     - Several random numbers generators are available through the bpp::RandomNumberFactory interface, and the
 *       bpp::RandomTools static class provides useful methods to deal with random numbers.
 *     .
 *   .
 * .
 */

namespace bpp
{

/**
 * @brief The Clonable interface (allow an object to be cloned).
 *
 * A clone is a deep (or hard) copy of an object.
 * This interface provides a method that dynamically creates a copy of itself
 * and send a pointer toward it.
 *
 * This method allows an object to be copied when you do not now its class.
 */
class Clonable
{
  public:
    
    Clonable() {}
  
    virtual ~Clonable() {}
  
  public:
    
    /**
     * @brief Create a copy of this object and send a pointer to it.
     *
     * @return A pointer toward the copy object.
     */
    virtual Clonable * clone() const = 0;
};

} //end of namespace bpp.

#endif  //_CLONABLE_H_

