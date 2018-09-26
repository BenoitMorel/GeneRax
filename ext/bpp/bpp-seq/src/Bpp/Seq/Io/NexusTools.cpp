//
// File: NexusTools.cpp
// Created by: Julien Dutheil
// Created on: Wed May 27 19:30 2009
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for sequences analysis.

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

#include "NexusTools.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Io/FileTools.h>

using namespace bpp;
using namespace std;

std::string NexusTools::getNextNonCommentLine(std::istream& input)
{
  string line = TextTools::removeSurroundingWhiteSpaces(FileTools::getNextLine(input));
  bool test = true;
  unsigned int countOpen = 0;
  unsigned int countClosed = 0;
  while(test)
  {
    if (line[0] == '[')
    {
      countOpen++;
    }
    if (line[line.size() - 1] == ']')
    {
      countClosed++;
    }
    if(countOpen > 0)
      line = TextTools::removeSurroundingWhiteSpaces(FileTools::getNextLine(input));
    if(countOpen == countClosed)
      test = false;
  }
  return line;
}


bool NexusTools::getNextCommand(std::istream& input, std::string& name, std::string& arguments, bool lineBrk)
noexcept(false)
{
	// Checking if the stream is readable
	if (! input) { throw IOException ("NexusTools::getNextCommand(). Failed to read from stream"); }
	
  string line = TextTools::removeSurroundingWhiteSpaces(getNextNonCommentLine(input));
  if (TextTools::startsWith(line, "BEGIN"))
  {
	  return false;
  }

  // Check if the command stands on one line:
  bool commandComplete = TextTools::endsWith(line, ";");
  if (commandComplete)
    line = line.substr(0, line.size() - 1);
  // Get the command name, as the first block:
  string::size_type limit = line.find(" ");
  if (limit == string::npos)
  {
    name = line;
    arguments = "";
    if (commandComplete)
    {
      //Command with no argument:
      return true;
    }
  }
  else
  {
    name = line.substr(0, limit);
    arguments = line.substr(limit + 1);
  }
  //Then parse the next lines:
  while(!commandComplete)
  {
	  if (input.eof()) { throw IOException ("NexusTools::getNextCommand(). Reached end of file before the end of the command could be found"); }
    line = TextTools::removeSurroundingWhiteSpaces(getNextNonCommentLine(input));
    commandComplete = TextTools::endsWith(line, ";");
    if (commandComplete)
      line = line.substr(0, line.size() - 1);
    if(lineBrk)
      arguments += "\n";
    arguments += line;
  }
  return true;
}

