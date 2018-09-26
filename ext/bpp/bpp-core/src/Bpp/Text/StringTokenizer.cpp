//
// File: StringTokenizer.cpp
// Author : Julien Dutheil
//          Sylvain Gaillard
// Last modification : Monday September 20 2004
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

#include "StringTokenizer.h"

using namespace bpp;
using namespace std;

StringTokenizer::StringTokenizer(const std::string& s, const std::string& delimiters, bool solid, bool allowEmptyTokens):
  tokens_(),
  splits_(),
  currentPosition_(0)
{
	if (!solid)
  {
    string::size_type index = s.find_first_not_of(delimiters, 0);
		while( index != s.npos)
    {
      string::size_type newIndex = s.find_first_of(delimiters, index);
			if (newIndex != s.npos)
      {
				tokens_.push_back(s.substr(index, newIndex - index));
				if (!allowEmptyTokens) index = s.find_first_not_of(delimiters, newIndex);
        else                   index = newIndex + 1;
        splits_.push_back(s.substr(newIndex, index - newIndex));
			}
      else
      {
				tokens_.push_back(s.substr(index));
				index = newIndex;
			}
		}
	}
	else
  {
    string::size_type index = 0;
		while (index != s.npos)
    {
      string::size_type newIndex = s.find(delimiters, index);
			if (newIndex != s.npos)
      {
				tokens_.push_back(s.substr(index, newIndex - index));
				if (!allowEmptyTokens)
        {
          index = newIndex + delimiters.size();
          while (index != string::npos && s.substr(index, delimiters.size()) == delimiters)
            index += delimiters.size();
        }
        else index = newIndex + delimiters.size();
				splits_.push_back(s.substr(newIndex, index - newIndex));
			}
      else
      {
				tokens_.push_back(s.substr(index));
				index = newIndex;
			}
		}
	}
}

void StringTokenizer::removeEmptyTokens()
{
  for (size_t i = tokens_.size(); i > currentPosition_; i--)
  {
    if (tokens_[i - 1] == "") tokens_.erase(tokens_.begin() + static_cast<ptrdiff_t>(i - 1));
  }
}

std::string StringTokenizer::unparseRemainingTokens() const
{
  string s;
  for (size_t i = currentPosition_; i < tokens_.size() - 1; ++i) {
    s += tokens_[i] + splits_[i];
  }
  if (numberOfRemainingTokens() > 0)
    s += tokens_.back();
  return s;
}

