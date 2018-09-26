//
// File: AttributesTools.h
// Created by: Julien Dutheil
// Created on: 2003
//

/*
   Copyright or © or Copr. Bio++ Development Tools, (November 17, 2004)

   This software is a computer program whose purpose is to provide basal and
   utilitary classes. This file belongs to the Bio++ Project.

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

#ifndef _ATTRIBUTES_TOOLS_H_
#define _ATTRIBUTES_TOOLS_H_

#include "../Exceptions.h"

// From the STL:
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

namespace bpp
{
/**
 * @brief Some functions to deal with attributes, i.e. parameters passed to a program/method.
 *
 * These methods allows you to retrieve attributes from command line arguments or
 * from a file.
 * The underlying syntax is <code> attributeName = argument </code>.
 * Here the delimiter char is '=', but another character may be used.
 *
 * In files, shell comments: <code> # my comment line here </code>,
 * C comments: <code> / * my comment block here * / </code> (but multiline not supported!)
 * and C++ comments: <code> // my comment line here </code> are allowed, and ignored while parsing.
 *
 * Lines may be broken, using the bash character '\' at the end of the line:
 * @code
 * optionfile=/home/foo/\
 *    bar.txt
 * @endcode
 * will be read as
 * @code
 * optionfile=/home/foo/bar.txt
 * @endcode
 *
 * Attributes are stored as a map<string, string>, with attributes names as keys,
 * and arguments as values.
 *
 * Here is an example of use.
 * This piece of code typically is at the begining of the main function.
 * It uses the FileTools and ApplicationTools classes, for checking file existence and displaying messages respectively.
 * @code
 * // Get the parameters from command line:
 * map<string, string> cmdParams = AttributesTools::getAttributesMap(
 *     AttributesTools::getVector(argc, argv), "=");
 *
 * // Look for a specified file with parameters:
 * int main(int argc, char *argv[]) {
 *   map<string, string> params;
 *   if(cmdParams.find("param") != cmdParams.end()) {
 *     string file = cmdParams["param"];
 *     if(!FileTools::fileExists(file)) {
 *       ApplicationTools::displayError("Parameter file not found.");
 *       exit(-1);
 *     } else {
 *       params = AttributesTools::getAttributesMapFromFile(file, "=");
 *       // Actualize attributes with the ones passed to command line:
 *       AttributesTools::actualizeAttributesMap(params, cmdParams);
 *     }
 *   } else {
 *     params = cmdParams;
 *   }
 * ...
 * @endcode
 * These pieces of code does the following:
 * - get all parameters from the command line and store them in a map,
 * - check if some parameter is called 'param' or 'params'. If so, look for the file
 *   given as value and try to read some parameters from it.
 * - If an parameter file was found, update the parameter in it with those from
 *   the command line. This implies that when a parameter is found in both the
 *   command line and the option file, the value from the command line will be
 *   retained.
 *
 * Support for variable is also available.
 * A variable character is specified (typically '$()') and may be used to de-reference any argument.
 * For instance,
 * @code
 * data=LSU
 * file=$(data).out
 * @endcode
 *
 * will be equivalent to
 * @code
 * data=LSU
 * file=LSU.out
 * @endcode
 */
class AttributesTools
{
public:
  AttributesTools() {}
  virtual ~AttributesTools() {}

  /**
   * @brief Get attributes a vector of strings from command line arguments.
   *
   * @param argc The number of arguments.
   * @param argv The array with all arguments.
   * @return A vector with all arguments as strings.
   */
  static std::vector<std::string> getVector(int argc, char* argv[]);

  /**
   * @brief Get an attribute map from a vector of arguments.
   *
   * This method also resolve all variable calls.
   *
   * @param argv      The vector of arguments.
   * @param delimiter The string that separates attribute names from arguments (for instance '=').
   * @return          The attribute map.
   */
  static std::map<std::string, std::string> getAttributesMap(
    const std::vector<std::string>& argv,
    const std::string& delimiter = "=");

  /**
   * @brief Get an attribute map from a vector of arguments.
   *
   * This method also resolve all variable calls.
   *
   * @param argv      The vector of arguments.
   * @param am        The attribute map to fill.
   * @param delimiter The string that separates attribute names from arguments (for instance '=').
   */
  static void getAttributesMap(
    const std::vector<std::string>& argv,
    std::map<std::string, std::string>& am,
    const std::string& delimiter = "=");

  /**
   * @brief Get an attribute map from a file.
   *
   * @param file      The file with all arguments.
   * @param delimiter The string that separates attribute names from arguments (for instance '=').
   * @return An attribute map.
   */
  static std::map<std::string, std::string> getAttributesMapFromFile(
    const std::string& file,
    const std::string& delimiter);

  /**
   * @brief Get an attribute map from a file.
   *
   * @param file      The file with all arguments.
   * @param params    An attribute map to fill.
   * @param delimiter The string that separates attribute names from arguments (for instance '=').
   */
  static void getAttributesMapFromFile(
    const std::string& file,
    std::map<std::string, std::string>& params,
    const std::string& delimiter);

  /**
   * @brief Actualizes an attribute map with another.
   *
   * All fields in map 2 will be added to map 1.
   * If two attributes have the same name, then map 1 value will be overwritten.
   *
   * @param attMap The attributes map.
   * @param atts   The attributes to add to the map.
   */
  static void actualizeAttributesMap(
    std::map<std::string, std::string>& attMap,
    const std::map<std::string, std::string>& atts);

  /**
   * @brief Resolve the variables.
   *
   * If used prior to the actualizeAttributesMap, this function will make the
   * variables 'local', whereas using them after will make them 'global'.
   *
   * @param am The attributes map.
   * @param varCode   The code that defines variable recalls.
   * @param varBeg    Variables begin name code.
   * @param varEnd    Variables end name code.
   * @throw Exception If there is a syntax error.
   */
  static void resolveVariables(std::map<std::string, std::string>& am,
                               char varCode = '$',
                               char varBeg = '(',
                               char varEnd = ')') noexcept(false);

  /**
   * @brief Global function that reads all parameters from command line and files,
   * and set the values in a map.
   *
   * @param args Number of arguments, as passed to the main function.
   * @param argv Array of values, as passed to the main function.
   * @return An attributes map.
   * @throw Exception in case an option file is not found.
   */
  static std::map<std::string, std::string> parseOptions(int args, char** argv)
  noexcept(false);

private:
  /**
   * @brief Remove comments from a string.
   *
   * @param s     The string to parse.
   * @param begin Comments front delimiter.
   * @param end   Comments end delimiter.
   */
  static std::string removeComments(const std::string& s, const std::string& begin, const std::string& end);
};
} // end of namespace bpp.

#endif // _ATTRIBUTES_TOOLS_H_
