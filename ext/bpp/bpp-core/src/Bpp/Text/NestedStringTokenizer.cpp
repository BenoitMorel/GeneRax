//
// File: NestedStringTokenizer.cpp
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

#include "NestedStringTokenizer.h"
#include "TextTools.h"

using namespace bpp;

//From the STL:
#include <iostream>

using namespace std;

NestedStringTokenizer::NestedStringTokenizer(const std::string& s, const std::string& open, const std::string& end, const std::string& delimiters, bool solid)
noexcept(false):
  StringTokenizer()
{
  int blocks = 0;
  string cache = "";
  if (!solid)
  {
    string::size_type index = s.find_first_not_of(delimiters, 0);
    while (index != s.npos)
    {
      string::size_type newIndex = s.find_first_of(delimiters, index);
      bool endBlockFound = false;
      while (!endBlockFound)
      {
        if (newIndex != s.npos)
        {
          string token = s.substr(index, newIndex - index);
          blocks += static_cast<int>(TextTools::count(token, open)) - static_cast<int>(TextTools::count(token, end));
        
          if (blocks == 0)
          {
            tokens_.push_back(cache + token);
            cache = ""; //reset cache.
            index = s.find_first_not_of(delimiters, newIndex);
            endBlockFound = true;
          }
          else
          {
            // Ignore this token untill closing block is found
            cache += s.substr(index, newIndex - index + 1);
            index = newIndex + 1;
            newIndex = s.find_first_of(delimiters, index);
          }
        }
        else
        {
          string token = s.substr(index);
          blocks += static_cast<int>(TextTools::count(token, open)) - static_cast<int>(TextTools::count(token, end));
          if (blocks == 0)
          {
            tokens_.push_back(cache + token);
            cache = ""; //reset cache.
            index = newIndex;
            endBlockFound = true;
          }
          else throw Exception("NestedStringTokenizer (constructor). Unclosed block.");
        }
      }
    }
  }
  else
  {
    string::size_type index = 0;
    while (index != s.npos)
    {
      string::size_type newIndex = s.find(delimiters, index);
      bool endBlockFound = false;
      while (!endBlockFound)
      {
        if (newIndex != s.npos)
        {
          string token = s.substr(index, newIndex - index);
          blocks += static_cast<int>(TextTools::count(token, open)) - static_cast<int>(TextTools::count(token, end));
				  
          if (blocks == 0)
          {
            tokens_.push_back(cache + token);
            cache = ""; //reset cache.
            index = newIndex + delimiters.size();
            endBlockFound = true;
          }
          else
          {
            // Ignore this token untill closing block is found
            cache += s.substr(index, newIndex - index + 1);
            index = newIndex + 1;
            newIndex = s.find(delimiters, index);
          }
        }
        else
        {
          string token = s.substr(index);
          blocks += static_cast<int>(TextTools::count(token, open)) - static_cast<int>(TextTools::count(token, end));
          if (blocks == 0)
          {
            tokens_.push_back(cache + token);
            cache = ""; //reset cache.
            index = newIndex;
            endBlockFound = true;
          }
          else throw Exception("Unclosed block."); 
        }
      }
    }
  }
}

const std::string& NestedStringTokenizer::nextToken() noexcept(false)
{
  if (!hasMoreToken()) throw Exception("No more token in nested tokenizer.");
  return tokens_[currentPosition_++];
}

