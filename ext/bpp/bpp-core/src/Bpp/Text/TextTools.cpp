//
// File: TextTools.cpp
// Created by: Julien Dutheil
// Created on: Fri Aug  8 12:57:50 2003
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

#include "TextTools.h"
#include "../Exceptions.h"

using namespace bpp;

#include <ctype.h>
#include <sstream>
#include <iomanip>

using namespace std;

/******************************************************************************/

bool TextTools::isEmpty(const std::string& s)
{
  for (unsigned int i = 0; i < s.size(); i++)
  {
    char c = s[i];
    if (c != ' ' && c != '\n' && c != '\t')
      return false;
  }
  return true;
}

/******************************************************************************/

std::string TextTools::toUpper(const std::string& s)
{
  string result = "";
  for (size_t i = 0; i < s.size(); i++)
  {
    result += static_cast<char>(toupper(static_cast<int>(s[i])));
  }
  return result;
}

/******************************************************************************/

std::string TextTools::toLower(const std::string& s)
{
  string result = "";
  for (size_t i = 0; i < s.size(); i++)
  {
    result += static_cast<char>(tolower(static_cast<int>(s[i])));
  }
  return result;
}

/******************************************************************************/

bool TextTools::isWhiteSpaceCharacter(char c)
{
  return (c == ' ')
         || (c == '\t')
         || (c == '\n')
         || (c == '\r')
         || (c == '\f');
}

/******************************************************************************/

std::string TextTools::removeWhiteSpaces(const std::string& s)
{
  // Copy sequence
  string st (s);

  // For all sequence's characters
  for (unsigned int i = 0; i < st.size(); i++)
  {
    if (isWhiteSpaceCharacter(st[i]))
    {
      st.erase(st.begin() + i); // Remove character
      i--;
    }
  }

  // Send result
  return st;
}

/******************************************************************************/

std::string TextTools::removeFirstWhiteSpaces(const std::string& s)
{
  // Copy sequence
  string st (s);

  while (st.size() > 0 && isWhiteSpaceCharacter(st[0]))
  {
    st.erase(st.begin());
  }

  // Send result
  return st;
}

/******************************************************************************/

std::string TextTools::removeLastWhiteSpaces(const std::string& s)
{
  // Copy sequence
  string st (s);

  while (st.size() > 0 && isWhiteSpaceCharacter(st[st.size() - 1]))
  {
    st.erase(st.end() - 1);
  }

  // Send result
  return st;
}

/******************************************************************************/

std::string TextTools::removeSurroundingWhiteSpaces(const std::string& s)
{
  return removeFirstWhiteSpaces(removeLastWhiteSpaces(s));
}

/******************************************************************************/

bool TextTools::isNewLineCharacter(char c)
{
  return (c == '\n')
         || (c == '\r');
}

/******************************************************************************/

std::string TextTools::removeNewLines(const std::string& s)
{
  // Copy string
  string st (s);

  // For all string's characters
  for (unsigned int i = 0; i < st.size(); i++)
  {
    if (isNewLineCharacter(st[i]))
    {
      st.erase(st.begin() + i); // Remove character
      i--;
    }
  }

  // Send result
  return st;
}

/******************************************************************************/

std::string TextTools::removeLastNewLines(const std::string& s)
{
  // Copy string
  string st (s);

  while (st.size() > 0 && isNewLineCharacter(st[st.size() - 1]))
  {
    st.erase(st.end() - 1);
  }

  // Send result
  return st;
}

/******************************************************************************/

bool TextTools::isDecimalNumber(char c)
{
  if (c == '0' || c == '1' || c == '2' || c == '3' || c == '4'
      || c == '5' || c == '6' || c == '7' || c == '8' || c == '9')
    return true;
  else
    return false;
}

/******************************************************************************/

bool TextTools::isDecimalNumber(const std::string& s, char dec, char scientificNotation)
{
  if (isEmpty(s))
    return false;
              
  size_t sepCount = 0;
  size_t sciCount = 0;
  size_t i = 0;
  if (s[0] == '-') i = 1;
  for (; i < s.size(); ++i)
  {
    char c = s[i];
    if (c == dec)
      sepCount++;
    else if (c == scientificNotation) {
      sciCount++;
      if (i == s.size() - 1) return false; //Must be sthg after scientific notation.
      c = s[i + 1];
      if (c == '-' || c == '+') i++;
      if (i == s.size() - 1) return false; //Must be sthg after scientific notation.
      if (sepCount == 0) sepCount = 1; //We do not want any dec in the exponent.
    } else if (!isDecimalNumber(c))
      return false;
    if (sepCount > 1 || sciCount > 1)
      return false;
  }
  return true;
}

/******************************************************************************/

bool TextTools::isDecimalInteger(const std::string& s, char scientificNotation)
{
  if (isEmpty(s))
    return false;
              
  size_t sciCount = 0;
  size_t i = 0;
  if (s[0] == '-') i = 1;
  for (; i < s.size(); ++i)
  {
    char c = s[i];
    if (c == scientificNotation) {
      sciCount++;
      if (i == s.size() - 1) return false; //Must be sthg after scientific notation.
      c = s[i + 1];
      if (c == '-') return false; //Not an integer then!
      if (c == '+') i++;
      if (i == s.size() - 1) return false; //Must be sthg after scientific notation.
    } else if (!isDecimalNumber(c))
      return false;
    if (sciCount > 1)
      return false;
  }
  return true;
}

/******************************************************************************/

std::string TextTools::toString(int i)
{
  ostringstream oss;
  oss << i;
  return oss.str();
}

/******************************************************************************/

std::string TextTools::toString(char c)
{
  ostringstream oss;
  oss << c;
  return oss.str();
}

/******************************************************************************/

std::string TextTools::toString(double d, int precision)
{
  ostringstream oss;
  oss << setprecision(precision) << d;
  return oss.str();
}

/******************************************************************************/

int TextTools::toInt(const std::string& s, char scientificNotation)
{
  if (!isDecimalInteger(s, scientificNotation)) throw Exception("TextTools::toInt(). Invalid number specification: " + s);
  istringstream iss(s);
  int i;
  iss >> i;
  return i;
}

/******************************************************************************/

double TextTools::toDouble(const std::string& s, char dec, char scientificNotation)
{
  if (!isDecimalNumber(s, dec, scientificNotation)) throw Exception("TextTools::toDouble(). Invalid number specification: " + s);
  istringstream iss(s);
  double d;
  iss >> d;
  return d;
}

/******************************************************************************/

std::string TextTools::resizeRight(const std::string& s, size_t newSize, char fill)
{
  if (s.size() > newSize)
    return s.substr(0, newSize);
  else
    return s + std::string(newSize - s.size(), fill);
}

/******************************************************************************/

std::string TextTools::resizeLeft(const std::string& s, size_t newSize, char fill)
{
  if (s.size() > newSize)
    return s.substr(s.size() - newSize);
  else
    return string(newSize - s.size(), fill) + s;
}

/******************************************************************************/

std::vector<std::string> TextTools::split(const std::string& s, size_t n)
{
  vector<string> v;
  string tmp = s;
  while (tmp.size() > n)
  {
    v.push_back(tmp.substr(0, n));
    tmp = tmp.substr(n);
  }
  v.push_back(tmp);
  return v;
}

/******************************************************************************/

std::string TextTools::removeSubstrings(const std::string& s, char blockBeginning, char blockEnding)
{
  string t = "";
  int blockCount = 0;
  size_t begPos = 0;
  for (size_t i = 0; i < s.size(); i++)
  {
    char current = s[i];
    if (current == blockBeginning)
    {
      blockCount++;
      t += s.substr(begPos, i - begPos);
    }
    else if (current == blockEnding)
    {
      blockCount--;
      if (blockCount == 0)
      {
        begPos = i + 1;
      }
      else if (blockCount < 0)
        throw Exception("TextTools::removeSubstrings(). " +
                        string("Ending block character without corresponding beginning one at position ") + toString((int)i) + ".");
    }
  }
  t += s.substr(begPos);
  return t;
}

/******************************************************************************/

std::string TextTools::removeSubstrings(const std::string& s, char blockBeginning, char blockEnding, std::vector<string>& exceptionsBeginning, std::vector<string>& exceptionsEnding)
{
  string t = "";
  int blockCount = 0;
  size_t begPos = 0;
  for (size_t i = 0; i < s.size(); i++)
  {
    char current = s[i];
    if (current == blockBeginning)
    {
      bool except = false;
      for (size_t j = 0; j < exceptionsBeginning.size(); j++)
      {
        size_t pos = exceptionsBeginning[j].find(blockBeginning);
        if (pos != string::npos) {
          size_t left = i - pos;
          size_t right = i + exceptionsBeginning[j].length() - pos;
          if ((right < s.length() - 1) && (hasSubstring (s.substr(left, right), exceptionsBeginning[j])))
          {
            except = true;
            break;
          }
        }
      }
      if (!except)
      {
        blockCount++;
        t += s.substr(begPos, i - begPos);
      }
    }
    else if ( (current == blockEnding) && (blockCount > 0) )
    {
      for (size_t j = 0; j < exceptionsEnding.size(); j++)
      {
        size_t pos = exceptionsEnding[j].find(blockEnding);
        if (pos != string::npos) {
          size_t left = i - pos;
          size_t right = i + exceptionsEnding[j].length() - pos;
          if ((right < s.length() - 1 ) && (hasSubstring (s.substr(left, right), exceptionsEnding[j])))
          {
            break;
          }
        }
      }
      blockCount--;
      if (blockCount == 0)
      {
        begPos = i + 1;
      }
      else if (blockCount < 0)
        throw Exception("TextTools::removeSubstrings(). " +
                        string("Ending block character without corresponding beginning one at position ") + toString((int)i) + ".");
    }
  }
  t += s.substr(begPos);
  return t;
}

/******************************************************************************/

std::string TextTools::removeChar(const std::string& s, char c)
{
  // Copy sequence
  string st(s);

  // For all sequence's characters
  for (unsigned int i = 0; i < st.size(); i++)
  {
    if (st[i] == c)
    {
      st.erase(st.begin() + i); // Remove character
      i--;
    }
  }

  // Send result
  return st;
}

/******************************************************************************/

unsigned int TextTools::count(const std::string& s, const std::string& pattern)
{
  unsigned int count = 0;
  string::size_type index = s.find(pattern);
  while (index != string::npos)
  {
    count++;
    index = s.find(pattern, index + 1);
  }
  return count;
}

/******************************************************************************/

bool TextTools::endsWith(const std::string& s, const std::string& pattern)
{
  if (s.size() < pattern.size())
    return false;
  return s.substr(s.size() - pattern.size()) == pattern;
}

/******************************************************************************/

bool TextTools::hasSubstring(const std::string& s, const std::string& pattern)
{
  if (s.size() < pattern.size())
    return false;
  for (size_t i = 0; i < s.size() - pattern.size() + 1; ++i)
  {
    if (s.substr(i, pattern.size()) == pattern)
      return true;
  }
  return false;
}

/******************************************************************************/

void TextTools::replaceAll(std::string& target, const std::string& query, const std::string& replacement)
{
  if (query.empty())
    return;
  size_t pos = target.find(query);
  while (pos != string::npos) {
    target.replace(pos, query.length(), replacement);
    pos += replacement.length(); //We prevent recursivity! 
    pos = target.find(query, pos);
  }
}

/******************************************************************************/

