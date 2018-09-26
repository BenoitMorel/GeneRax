//
// File: KeyvalTools.cpp
// Created by: Julien Dutheil
// Created on: Mon May 13:16 CET 2009
//

/*
Copyright or © or Copr. Bio++ Development Team, (2009)

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

#include "KeyvalTools.h"
#include "NestedStringTokenizer.h"

//From the STL:
#include <memory>

using namespace bpp;
using namespace std;

void KeyvalTools::singleKeyval(const std::string& desc, std::string& key, std::string& val, const std::string& split) noexcept(false)
{
  string::size_type i = desc.find(split);
  if (i == string::npos)
    throw KeyvalException("Bad syntax! keyval should be of the form 'key" + split + "=value', found '" + desc + "'.");
  key = desc.substr(0, i);
  val = desc.substr(i+1);
}

void KeyvalTools::multipleKeyvals(const std::string& desc, std::map<std::string,std::string>& keyvals, const std::string& split, bool nested) noexcept(false)
{
  unique_ptr<StringTokenizer> st;
  if (nested)
    st.reset(new NestedStringTokenizer(desc, "(", ")", split));
  else
    st.reset(new StringTokenizer(desc, split));
  string key, val;
  vector<string> tokens;
  //Check tokens:
  string token;
  while (st->hasMoreToken())
  {
    token = st->nextToken();
    if (token == "=")
    {
      //We need to merge the next token with the last one:
      if (tokens.size() == 0)
        throw KeyvalException("Invalid syntax, found '=' without argument name.");
      if (!st->hasMoreToken())
        throw KeyvalException("Invalid syntax, found '=' without argument value.");
      string nextToken = st->nextToken();
      if (nextToken == "=")
        throw KeyvalException("Invalid syntax, found a double '='.");
      tokens[tokens.size() - 1] += "=" + nextToken;
    }
    else
    {
      tokens.push_back(token);
    }
  }
  for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); it++)
  {
    singleKeyval(*it, key, val);
    key=TextTools::removeSurroundingWhiteSpaces(key);
    val=TextTools::removeSurroundingWhiteSpaces(val);
    keyvals[key] = val;
  }
}

std::string KeyvalTools::changeKeyvals(const std::string& desc, const std::map<std::string, std::string>& newkeyvals, const std::string& split, bool nested) noexcept(false)
{
  string::size_type begin = desc.find_first_of("(");
  string::size_type end = desc.find_last_of(")");

  if (begin == string::npos && end == string::npos)
  {
    //Empty procedure:
    return desc;
  }
  if (begin == string::npos && end != string::npos)
    throw KeyvalException("Bad keyval procedure, missing opening parenthesis.");
  if (begin == string::npos && end != string::npos)
    throw KeyvalException("Bad keyval procedure, missing closing parenthesis.");
  
  if (!TextTools::isEmpty(desc.substr(end + 1)))
    throw KeyvalException("Bad keyval procedure, extra characters after closing parenthesis: " + desc.substr(end + 1));
  //Get the procedure name (without leading spaces):

  string newDesc= TextTools::removeFirstWhiteSpaces(desc.substr(0, begin))+"(";

  string desckv=desc.substr(begin + 1, end - begin - 1);

  unique_ptr<StringTokenizer> st;
  if (nested)
    st.reset(new NestedStringTokenizer(desckv, "(", ")", split));
  else
    st.reset(new StringTokenizer(desckv, split));
  string key, val;
  vector<string> tokens;
  //Check tokens:
  string token;
  
  while (st->hasMoreToken())
  {
    token = st->nextToken();
    if (token == "=")
    {
      //We need to merge the next token with the last one:
      if (tokens.size() == 0)
        throw KeyvalException("Invalid syntax, found '=' without argument name.");
      if (!st->hasMoreToken())
        throw KeyvalException("Invalid syntax, found '=' without argument value.");
      string nextToken = st->nextToken();
      if (nextToken == "=")
        throw KeyvalException("Invalid syntax, found a double '='.");
      tokens[tokens.size() - 1] += "=" + nextToken;
    }
    else
    {
      tokens.push_back(token);
    }
  }

  for (vector<string>::iterator it = tokens.begin(); it != tokens.end(); it++)
  {
    singleKeyval(*it, key, val);
    key=TextTools::removeSurroundingWhiteSpaces(key);
    if (it!=tokens.begin())
      newDesc+=split;

    map<string, string>::const_iterator iter=newkeyvals.find(key);
    
    if (iter!=newkeyvals.end())
      newDesc+=key+"="+iter->second;
    else
      newDesc+=*it;
  }

  newDesc+=")";
  
  return newDesc;
}

void KeyvalTools::parseProcedure(const std::string& desc, std::string& name, std::map<std::string, std::string>& args) noexcept(false)
{
  string::size_type begin = desc.find_first_of("(");
  string::size_type end = desc.find_last_of(")");
  
  if (begin == string::npos && end == string::npos)
  {
    //Empty procedure:
    name = desc;
    return;
  }
  if (begin == string::npos && end != string::npos)
    throw KeyvalException("Bad keyval procedure, missing opening parenthesis.");
  if (begin == string::npos && end != string::npos)
    throw KeyvalException("Bad keyval procedure, missing closing parenthesis.");
  
  if (!TextTools::isEmpty(desc.substr(end + 1)))
    throw KeyvalException("Bad keyval procedure, extra characters after closing parenthesis: " + desc.substr(end + 1));
  //Get the procedure name (without leading spaces):
  name = TextTools::removeFirstWhiteSpaces(desc.substr(0, begin));
  multipleKeyvals(desc.substr(begin + 1, end - begin - 1), args);
}

