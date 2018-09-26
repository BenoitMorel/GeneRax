/*
Copyright or © or Copr. Bio++ Development Team, (September 19, 2014)

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

#include "PhyloTree.h"

using namespace bpp;
using namespace std;

PhyloTree::PhyloTree(bool rooted) :
  AssociationTreeGlobalGraphObserver<PhyloNode,PhyloBranch>(rooted),
  name_("")
{
}

/*
PhyloTree::PhyloTree(const ParametrizablePhyloTree& tree) :
  AssociationTreeGlobalGraphObserver<PhyloNode,PhyloBranch>(tree),
  name_("")
{
}
*/


void PhyloTree::resetNodesId()
{
  std::vector<shared_ptr<PhyloNode> > nodes = getAllNodes();
  for (unsigned int i = 0; i < nodes.size(); i++)
  {
    setNodeIndex(nodes[i], i);
    
    if (hasFather(nodes[i]))
      setEdgeIndex(getEdgeToFather(nodes[i]),i);
  }
}

std::vector<std::string> PhyloTree::getAllLeavesNames() const
{
  vector<string> vn;
  
  vector<shared_ptr<PhyloNode> > vpn=getAllLeaves();

  for (vector<shared_ptr<PhyloNode> >::const_iterator it=vpn.begin(); it != vpn.end(); it++)
    vn.push_back((*it)->getName());
  
  return vn;
}

void PhyloTree::scaleTree(shared_ptr<PhyloNode> node, double factor)
{
  vector<shared_ptr<PhyloBranch> > branches = getSubtreeEdges(node);
  for (vector<shared_ptr<PhyloBranch> >::iterator currBranch = branches.begin(); currBranch != branches.end(); currBranch++)
  {
    if((*currBranch)->hasLength())
      (*currBranch)->setLength((*currBranch)->getLength()*factor);
    else
      throw PhyloBranchPException("PhyloTree::scaleTree : undefined length", (*currBranch).get());
  }
}
      
void PhyloTree::scaleTree(double factor)
{
  scaleTree(getRoot(), factor);
}

void PhyloTree::setBranchLengths(double l)
{
  vector<shared_ptr<PhyloBranch> > vpn=getAllEdges();

  for (vector<shared_ptr<PhyloBranch> >::const_iterator it=vpn.begin(); it != vpn.end(); it++)
    (*it)->setLength(l);
}

