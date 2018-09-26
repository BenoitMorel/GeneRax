//
// File: PhyloTreeExceptions.cpp
// Created by: Julien Dutheil
// Created on: Mon Nov  3 17:04:46 2003
//

/*
   Copyright or © or Copr. CNRS, (November 16, 2004)

   Julien.Dutheil@univ-montp2.fr

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

#include "PhyloTreeExceptions.h"
#include "PhyloNode.h"
#include "PhyloBranch.h"
#include "PhyloTree.h"

#include <Bpp/Text/TextTools.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

PhyloNodePException::PhyloNodePException(const std::string& text, const PhyloTree& tree, const shared_ptr<PhyloNode> node) :
  PhyloNodeException(text, tree.getNodeIndex(node)), node_(node.get())
{}

/******************************************************************************/

PhyloNodePException::PhyloNodePException(const std::string& text, const PhyloNode* node) :
  PhyloNodeException(text, 0), node_(node)
{}

/******************************************************************************/

PhyloNodeNotFoundException::PhyloNodeNotFoundException(const std::string& text, const std::string& id) :
  Exception("NodeNotFoundException: " + text + "(" + id + ")"),
  id_(id) {}

PhyloNodeNotFoundException::PhyloNodeNotFoundException(const std::string& text, int id) :
  Exception("NodeNotFoundException: " + text + "(" + TextTools::toString(id) + ")"),
  id_(TextTools::toString(id)) {}

/******************************************************************************/

PhyloBranchPException::PhyloBranchPException(const std::string& text, const PhyloTree& tree, const shared_ptr<PhyloBranch> branch) :
  PhyloBranchException(text, tree.getEdgeIndex(branch)), branch_(branch.get())
{}

/******************************************************************************/

PhyloBranchPException::PhyloBranchPException(const std::string& text, const PhyloBranch* branch) :
  PhyloBranchException(text, 0), branch_(branch)
{}

/******************************************************************************/

PhyloBranchNotFoundException::PhyloBranchNotFoundException(const std::string& text, const std::string& id) :
  Exception("BranchNotFoundException: " + text + "(" + id + ")"),
  id_(id) {}

/******************************************************************************/

PhyloBranchNotFoundException::PhyloBranchNotFoundException(const std::string& text, int id) :
  Exception("BranchNotFoundException: " + text + "(" + TextTools::toString(id) + ")"),
  id_(TextTools::toString(id)) {}

/******************************************************************************/

PhyloTreeException::PhyloTreeException(const std::string& text, const PhyloTree* tree) :
  Exception("PhyloTreeException: " + text + (tree != 0 ? "(" + tree->getName() + ")" : "")),
  tree_(tree) {}

/******************************************************************************/

UnrootedPhyloTreeException::UnrootedPhyloTreeException(const std::string& text, const PhyloTree* tree) :
  PhyloTreeException("UnrootedPhyloTreeException: " + text, tree) {}

/******************************************************************************/

