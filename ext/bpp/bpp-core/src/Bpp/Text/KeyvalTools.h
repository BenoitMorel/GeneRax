//
// File: KeyvalTools.h
// Created by: Julien Dutheil
// Created on: Mon May 13:16 CET 2009
//

/*
Copyright or © or Copr. Bio++ Development Tools, (2009)

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

#ifndef _KEYVALTOOLS_H_
#define _KEYVALTOOLS_H_

#include "StringTokenizer.h"
#include "TextTools.h"
#include "../Exceptions.h"

//From the STL:
#include <map>

namespace bpp
{

/**
 * @brief Exception thrown by the Keyval parser.
 */
class KeyvalException :
  public Exception
{
  public:
    KeyvalException(const std::string& message) : Exception(message) {}
};

/**
 * @brief Tools to deal with the keyval syntax.
 *
 * This class contains method to deal with parameter=value syntax procedure.
 * A keyval procedure takes the form 
 * @code
 * proc(p1=v1,p2=v2,p3=v3,etc)
 * @endcode
 * where 'p' are parameter names, and 'v' are the corresponding values.
 * These values can be nested keyval procedures.
 */
class KeyvalTools
{
  public:
    KeyvalTools();
    virtual ~KeyvalTools();

  public:
    /**
     * @brief Split a string into a key and a value (General purpose function).
     *
     * @param desc  [in]  A string descibing the keyval, with format key=val (space are considered normal character, that's up to you to deal with that afterward!).
     * @param key   [out] Will contain the text of the key.
     * @param val   [out] Will contain the text of the value.
     * @param split [in]  The delimiter. Default is '=' but ':' can be used.
     * @throw KeyvalException If the syntax describing the keyval is not correct.
     */

  static void singleKeyval(const std::string& desc, std::string& key, std::string& val, const std::string& split = "=") noexcept(false);
    
    /**
     * @brief Split a string into several keys and corresponding values (General purpose function).
     *
     * @param desc [in]  A string descibing the keyval, with format key1=val1,key2=val2,etc (space are considered normal character, that's up to you to deal with that afterward!).
     * @param keyvals [out] Will contain the text of the keys and their corresponding values.
     * @param split [in] The keyval delimiter.the default is a coma, but a space character can be used for instance.
     * @param nested [in] Tell if nested keyval procedures are expected.
     * @throw KeyvalException If the syntax describing the keyval is not correct.
     */
  
    static void multipleKeyvals(const std::string& desc, std::map<std::string, std::string>& keyvals, const std::string& split = ",", bool nested = true) noexcept(false);

  /**
   * @brief Change several keys to new corresponding values (General
   * purpose function).
   *
   * @param desc [in]  A string descibing the keyval, with format
   * key1=val1,key2=val2,etc (space are considered normal character,
   * that's up to you to deal with that afterwards!).
   * @param newkeyvals [in] contains the text of the keys to be changed
   * and their new corresponding values. If a key is not in desc, it
   * is not added.
   * @param split [in] The keyval delimiter. The default is a coma,
   * but a space character can be used for instance.
   * @param nested [in] Tell if nested keyval procedures are expected.
   * @return the string with the changed values.
   * @throw KeyvalException If the syntax describing the keyval is not correct.
   */
  
  static std::string changeKeyvals(const std::string& desc, const std::map<std::string, std::string>& newkeyvals, const std::string& split = ",", bool nested = true) noexcept(false);

    /**
     * @brief Parse (not recursively) a procedure string.
     *
     * @param desc [in]  A string descibing the keyval procedure.
     * @param name [out] Outputs the name of the procedure.
     * @param args [out] Fills a map with all keys and values for parameters.
     * @throw KeyvalException If the description is invalid (one parenthesis is missing for instance).
     */
    static void parseProcedure(const std::string& desc, std::string& name, std::map<std::string, std::string>& args) noexcept(false);
};

} //End of namespace bpp.

#endif  //_KEYVALTOOLS_H_

