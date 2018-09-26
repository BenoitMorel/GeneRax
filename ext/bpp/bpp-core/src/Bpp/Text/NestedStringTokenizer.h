//
// File: NestedStringTokenizer.h
// Author : Julien Dutheil
// Last modification : Monday May 22 10:57 2006
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

  This software is a computer program whose purpose is to map data onto
  a sequence or a phylogenetic tree.

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

#ifndef _NESTEDSTRINGTOKENIZER_H_
#define _NESTEDSTRINGTOKENIZER_H_

//From the STL:
#include <deque>
#include <string>

#include "StringTokenizer.h"
#include "../Exceptions.h"

namespace bpp
{

  /**
   * @brief An improved tokenizer for strings.
   *
   * Splits a string according to a given (set of) delimiter(s).
   * Delimiters in certains blocks ({}, [], etc) are ignored.
   */
  class NestedStringTokenizer:
    public StringTokenizer
  {
  public:
		
    /**
     * @brief Build a new StringTokenizer from a string.
     *
     * @param s          The string to parse.
     * @param open       Opening block.
     * @param end        Ending block.
     * @param delimiters Chars that must be considered as delimiters.
     * @param solid      If true, delimiters is considered as a single bloc delimiter.
     */
    NestedStringTokenizer(const std::string& s, const std::string& open, const std::string& end, const std::string& delimiters = " \t\n\f\r", bool solid = false) noexcept(false);
		
    virtual ~NestedStringTokenizer() {}
	
  public:
		
    /**
     * @brief Get the next available token.
     * If no token is availbale, throw an Exception.
     *
     * @return The next token if there is one.
     */
    const std::string& nextToken() noexcept(false);


    /**
     * @brief This function is not supported for nested tokenizers.
     *
     * @return An empty string.
     */
    std::string unparseRemainingTokens() const { return ""; }
  };

} //end of namespace bpp;

#endif	//_NESTEDSTRINGTOKENIZER_H_

