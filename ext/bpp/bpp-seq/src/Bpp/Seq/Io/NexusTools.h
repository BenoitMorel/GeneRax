//
// File: NexusTools.h
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

#ifndef _NEXUSTOOLS_H_
#define _NEXUSTOOLS_H_

// From the STL:
#include <iostream>
#include <Bpp/Exceptions.h>

namespace bpp
{

/**
 * @brief Tools for parsing Nexus files.
 *
 * The Nexus format is described in the following paper:
 * Maddison D, Swofford D, and Maddison W (1997), _Syst Biol_ 46(4):590-621
 *
 * @author Julien Dutheil
 */
class NexusTools
{
  public:
    /**
     * @param input The input stream.
     * @return A string containing the next line in the file wichi is not empty and is no a comment line.
     */
    static std::string getNextNonCommentLine(std::istream& input);

    /**
     * @brief parse the next command name within a block.
     *
     * @param input     [in]  The input stream.
     * @param name      [out] Will contain the name of the command. 
     * @param arguments [out] Will contain the arguments of the commans, as raw data. The arguments will not be parsed.
     * @param lineBrk   [in]  Tell is the line break should be preserved in the arguments.
     * @return Whether a command was found in the current block.
     * @throw IOException In case of bad format.
     */
    static bool getNextCommand(std::istream& input, std::string& name, std::string& arguments, bool lineBrk = true)
    noexcept(false);

};

} //end of namespace bpp.

#endif  //_NEXUSTOOLS_H_

