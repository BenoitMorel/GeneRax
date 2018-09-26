//
// File: NexusIOTree.cpp
// Created by: Julien Dutheil
// Created on: Wed May 27 19:06 2009
//

/*
  Copyright or © or Copr. CNRS, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

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

#include "NexusIoTree.h"
#include "Newick.h"
#include "../Tree/PhyloTree.h"

#include <Bpp/Text/TextTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/NestedStringTokenizer.h>
#include <Bpp/Numeric/VectorTools.h>

//From bpp-seq:
#include <Bpp/Seq/Io/NexusTools.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>
#include <sstream>

using namespace std;

/******************************************************************************/

const string NexusIOTree::getFormatName() const { return "Nexus"; }

/******************************************************************************/

const string NexusIOTree::getFormatDescription() const
{
  return string("Nexus format (trees only). ");
}

/******************************************************************************/
/* INPUT */
/******************************************************************************/

/******************************************************************************/

PhyloTree * NexusIOTree::readP(istream &in) const
{
  vector<PhyloTree*> trees;
  read(in, trees);
  if (trees.size() == 0)
    throw IOException("NexusIOTree::read(). No tree found in file.");
  for (size_t i = trees.size() - 1; i > 0; i--)
    delete trees[i];
  return trees[0];
}

/******************************************************************************/

void NexusIOTree::read(std::istream& in, std::vector<PhyloTree*>& trees) const
{
  // Checking the existence of specified file
  if (! in) { throw IOException ("NexusIOTree::read(). Failed to read from stream"); }
	
  //Look for the TREES block:
  string line = "";
  while (TextTools::toUpper(line) != "BEGIN TREES;")
  {
    if (in.eof())
      throw Exception("NexusIOTree::read(). No trees block was found.");
    line = TextTools::removeSurroundingWhiteSpaces(FileTools::getNextLine(in));
  }
  
  string cmdName = "", cmdArgs = "";
  bool cmdFound = NexusTools::getNextCommand(in, cmdName, cmdArgs, false);
  if (! cmdFound)
    throw Exception("NexusIOTree::read(). Missing tree command.");
  cmdName = TextTools::toUpper(cmdName);

  //Look for the TRANSLATE command:
  map<string, string> translation;
  bool hasTranslation = false;
  if (cmdName == "TRANSLATE")
  {
    //Parse translation:
    StringTokenizer st(cmdArgs, ",");
    while (st.hasMoreToken())
    {
      string tok = TextTools::removeSurroundingWhiteSpaces(st.nextToken());
      NestedStringTokenizer nst(tok, "'", "'", " \t");
      if (nst.numberOfRemainingTokens() != 2)
        throw Exception("NexusIOTree::read(). Unvalid translation description.");
      string name = nst.nextToken();
      string tln  = nst.nextToken();
      translation[name] = tln;
    }
    hasTranslation = true;
    cmdFound = NexusTools::getNextCommand(in, cmdName, cmdArgs, false);
    if (! cmdFound)
      throw Exception("NexusIOTree::read(). Missing tree command.");
    else
      cmdName = TextTools::toUpper(cmdName);
  }

  //Now parse the trees:
  while (cmdFound && cmdName != "END")
  {
    if (cmdName != "TREE")
      throw Exception("NexusIOTree::read(). Unvalid command found: " + cmdName);
    string::size_type pos = cmdArgs.find("=");
    if (pos == string::npos)
      throw Exception("NexusIOTree::read(). unvalid format, should be tree-name=tree-description");
    string description = cmdArgs.substr(pos + 1);

    Newick treeReader;

    istringstream ss(description + ";");    
    PhyloTree* tree = treeReader.readP(ss);
    
    //Now translate leaf names if there is a translation:
    //(we assume that all trees share the same translation! ===> check!)
    if (hasTranslation)
    {
      vector<shared_ptr<PhyloNode> > leaves = tree->getAllLeaves();
      for (size_t i = 0; i < leaves.size(); i++)
      {
        string name = leaves[i]->getName();
        if (translation.find(name) == translation.end())
        {
          throw Exception("NexusIOTree::read(). No translation was given for this leaf: " + name);
        }
        leaves[i]->setName(translation[name]);
      }
    }
    trees.push_back(tree);
    cmdFound = NexusTools::getNextCommand(in, cmdName, cmdArgs, false);
    if (cmdFound) cmdName = TextTools::toUpper(cmdName);
  }
}

/******************************************************************************/
/*  OUTPUT */
/******************************************************************************/

/******************************************************************************/

void NexusIOTree::write_(const PhyloTree& tree, ostream& out) const
{
  vector<const PhyloTree*> trees;
  trees.push_back(&const_cast<PhyloTree&>(tree));
  write(trees, out);
}

/******************************************************************************/

void NexusIOTree::write_(const vector<const PhyloTree*>& trees, ostream& out) const
{
  // Checking the existence of specified file, and possibility to open
  // it in write mode
  if (! out) { throw IOException ("NexusIOTree::write: failed to write to stream"); }
  
  out << "#NEXUS" << endl;
  out << endl;
  out << "BEGIN TREES;" << endl;

  //First, we retrieve all leaf names from all trees:
  vector<string> names;
  for (size_t i = 0; i < trees.size(); i++)
  {
    names = VectorTools::vectorUnion(names, trees[i]->getAllLeavesNames());
  }
  //... and create a translation map:
  map<string, size_t> translation;
  size_t code = 0;
  for (size_t i = 0; i < names.size(); i++)
  {
    translation[names[i]] = code++;
  }

  //Second we translate all leaf names to their corresponding code:
  vector<PhyloTree*> translatedTrees;
  for (size_t i = 0; i < trees.size(); i++)
  {
    vector<shared_ptr<PhyloNode> > leaves = trees[i]->getAllLeaves();
    
    PhyloTree* tree = trees[i]->clone();
    
    for (size_t j = 0; j < leaves.size(); j++)
    {
      tree->getNode(trees[i]->getNodeIndex(leaves[j]))->setName(TextTools::toString(translation[leaves[j]->getName()]));
    }
    translatedTrees.push_back(tree);
  }

  //Third we print the translation command:
  out << "  TRANSLATE";
  size_t count = 0;
  for (map<string, size_t>::iterator it = translation.begin(); it != translation.end(); it++)
  {
    out << endl << "    " << it->second << "\t" << it->first;
    count++;
    if (count < translation.size())
      out << ",";
  }
  out << ";";

  Newick treeWriter;

  //Finally we print all tree descriptions:
  for (size_t i = 0; i < trees.size(); i++)
  {
    out << endl << "  TREE tree" << (i+1) << " = ";
    treeWriter.write(*translatedTrees[i], out);
  }
  out << "END;" << endl;
  
  //Clean trees:
  for (size_t i = 0; i < translatedTrees.size(); i++)
  {
    delete translatedTrees[i];
  }

}

/******************************************************************************/


