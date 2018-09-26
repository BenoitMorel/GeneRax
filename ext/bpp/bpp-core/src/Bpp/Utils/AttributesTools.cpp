//
// File: AttributesTools.cpp
// Created by: Julien Dutheil
// Created on: 2003
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

// From the STL:
#include <cstdlib>
#include <string>
#include <iostream>
#include <fstream>
#include <algorithm>

using namespace std;

#include "AttributesTools.h"
#include "../App/ApplicationTools.h"
#include "../Text/TextTools.h"
#include "../Io/FileTools.h"

using namespace bpp;

/******************************************************************************/

std::vector<std::string> AttributesTools::getVector(int argc, char* argv[])
{
  size_t n = static_cast<size_t>(argc);
  vector<string> result(n);
  for (size_t i = 1; i < n; ++i)
  {
    result[i] = string(argv[i]);
  }
  // Ignore first argc which is the program name!
  return result;
}

/******************************************************************************/

std::map<std::string, std::string> AttributesTools::getAttributesMap(
  const std::vector<std::string>& argv,
  const std::string& delimiter)
{
  map<string, string> am;
  getAttributesMap(argv, am, delimiter);
  return am;
}

/******************************************************************************/

void AttributesTools::getAttributesMap(
  const std::vector<std::string>& argv,
  std::map<std::string, std::string>& am,
  const std::string& delimiter)
{
  vector<string> argv2(argv.size());
  // First make a few cleaning:
  for (size_t i = 0; i < argv.size(); i++)
  {
    // Make a few corrections first:
    string arg = removeComments(argv[i], string("#"), string("\n")); // remove shell comments.
    arg = removeComments(arg, string("//"), string("\n")); // remove C simple comments.
    arg = removeComments(arg, string("/*"), string("*/")); // remove C multiple comments.
    arg = TextTools::removeWhiteSpaces(arg);
    argv2[i] = arg;
  }
  // Now parse arguments:
  for (size_t i = 0; i < argv.size(); i++)
  {
    string arg = argv2[i];
    if (arg == "")
      continue;  // Skipping void line.
    while (arg[arg.size() - 1] == '\\')
    {
      // Splitted line
      i++;
      arg = arg.substr(0, arg.length() - 1) + argv2[i];
    }
    // Parsing:
    string::size_type limit = arg.find(delimiter, 0);
    if (limit == string::npos)
    {
      // Invalid parameter
      (*ApplicationTools::warning << "WARNING!!! Parameter '" << arg << "' has been ignored.").endLine();
    }
    else
    {
      string name  = string(arg.begin(), arg.begin() + static_cast<ptrdiff_t>(limit));
      string value = string(arg.begin() + static_cast<ptrdiff_t>(limit + delimiter.size()), arg.end());
      // if ((name == "param") || (name == "params"))
      // {
      //   //Recursive inclusion:
      //   getAttributesMapFromFile(value, am, delimiter);
      // }
      // else
      am[name] = value;
    }
  }
}

/******************************************************************************/

void AttributesTools::getAttributesMapFromFile(
  const std::string& file,
  std::map<std::string, std::string>& params,
  const std::string& delimiter)
{
  cout << "Parsing file " << file << " for options." << endl;
  ifstream input(file.c_str(), ios::in);
  vector<string> lines = FileTools::putStreamIntoVectorOfStrings(input);
  getAttributesMap(lines, params, delimiter);
}

/******************************************************************************/

std::map<std::string, std::string> AttributesTools::getAttributesMapFromFile(
  const std::string& file,
  const std::string& delimiter)
{
  map<string, string> params;
  getAttributesMapFromFile(file, params, delimiter);
  return params;
}

/******************************************************************************/

void AttributesTools::actualizeAttributesMap(
  std::map<std::string, std::string>& attMap,
  const std::map<std::string, std::string>& atts)
{
  for (map<string, string>::const_iterator i = atts.begin(); i != atts.end(); i++)
  {
    if ((i->first != "param") && (i->first != "params"))
      attMap[i->first] = i->second;
  }
}

/******************************************************************************/

void AttributesTools::resolveVariables(
  std::map<std::string, std::string>& am,
  char varCode,
  char varBeg,
  char varEnd)
noexcept(false)
{
  // Now resolve any variable:
  for (map<string, string>::iterator it = am.begin(); it != am.end(); it++)
  {
    string value = it->second;
    string::size_type index1 = value.find(TextTools::toString(varCode) + TextTools::toString(varBeg));
    while (index1 != string::npos)
    {
      string::size_type index2 = value.find(TextTools::toString(varEnd), index1);
      if (index2 != string::npos)
      {
        string varName  = value.substr(index1 + 2, index2 - index1 - 2);
        map<string, string>::iterator varIt = am.find(varName);
        string varValue = "";
        if (varIt == am.end())
        {
          if (ApplicationTools::error)
            (*ApplicationTools::error << "Variable '" << varName << "' is undefined and was ignored.").endLine();
        }
        else
        {
          varValue = varIt->second;
        }
        // Modify original field:
        string newValue = value.substr(0, index1) + varValue + value.substr(index2 + 1);
        it->second = newValue;
      }
      else
        throw Exception("Syntax error, variable name is not closed.");
      value = it->second;
      index1 = value.find(TextTools::toString(varCode) + TextTools::toString(varBeg));
    }
  }
}

/******************************************************************************/

std::string AttributesTools::removeComments(
  const std::string& s,
  const std::string& begin,
  const std::string& end)
{
  string r = s;
  string::size_type last = 0;
  do
  {
    string::size_type first = r.find(begin, last);
    if (first == string::npos)
      return r;  // No shell comment.
    // else:
    last = r.find(end, first);
    if (last == string::npos)
    {
      r.erase(r.begin() + static_cast<ptrdiff_t>(first), r.end());
    }
    else
    {
      r.erase(r.begin() + static_cast<ptrdiff_t>(first), r.begin() + static_cast<ptrdiff_t>(last));
    }
  }
  while (last != string::npos);
  return r;
}

/******************************************************************************/

std::map<std::string, std::string> AttributesTools::parseOptions(int args, char** argv)
noexcept(false)
{
  // Get the parameters from command line:
  map<string, string> cmdParams = AttributesTools::getAttributesMap(
    AttributesTools::getVector(args, argv), "=");

  // Look for a specified file with parameters:
  map<string, string> params;
  std::map<std::string, std::string>::iterator it;

  if (cmdParams.find("param") != cmdParams.end())
  {
    string file = cmdParams["param"];
    if (!FileTools::fileExists(file))
    {
      throw Exception("AttributesTools::parseOptions(). Parameter file not found.");
    }
    else
    {
      params = getAttributesMapFromFile(file, "=");
      // Actualize attributes with ones passed to command line:
      actualizeAttributesMap(params, cmdParams);
    }
  }
  else
  {
    params = cmdParams;
  }
  // Resolve variables:
  resolveVariables(params);

  std::vector<string> mapfile;
  std::vector<string>::iterator imapfile;
  string file;

  while (true)
  {
    it = params.find("param");
    if (it != params.end())
    {
      file = it->second;
      if (std::find(mapfile.begin(), mapfile.end(), file) == mapfile.end())
      {
        params.erase(it);
        mapfile.push_back(file);
        getAttributesMapFromFile(file, params, "=");
        resolveVariables(params);
        continue;
      }
      else
        throw Exception("parsing error : Already used file " + file);
    }
    it = params.find("params");
    if (it != params.end())
    {
      file = it->second;
      if (find(mapfile.begin(), mapfile.end(), file) == mapfile.end())
      {
        params.erase(it);
        mapfile.push_back(file);
        getAttributesMapFromFile(file, params, "=");
        resolveVariables(params);
        continue;
      }
      else
        throw Exception("parsing error : Already used file " + file);
    }
    break;
  }

  return params;
}

/******************************************************************************/
