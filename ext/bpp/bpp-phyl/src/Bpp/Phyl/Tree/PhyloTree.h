//
// File: PhyloTree.h
// Created by: Laurent Guéguen
// Created on: dimanche 24 juillet 2016, à 19h 58
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

#ifndef _PHYLOTREE_H_
#define _PHYLOTREE_H_

#include <Bpp/Graph/AssociationTreeGraphImplObserver.h>

#include "PhyloNode.h"
#include "PhyloBranch.h"


namespace bpp
{
  /**
   * @brief Defines a Phylogenetic Tree based on a TreeGraph & its
   * associationObserver.
   *
   * @author Thomas Bigot
   */
  
  class ParametrizablePhyloTree;
    
  class PhyloTree:
    public AssociationTreeGlobalGraphObserver<PhyloNode,PhyloBranch>
  {
  private:
      std::string name_;
      
    public:

      PhyloTree(bool rooted = false);

      //PhyloTree(const ParametrizablePhyloTree& tree);

      PhyloTree* clone() const
      {
        return new PhyloTree(*this);
      }
      
      /**
       * @brief Tree name.
       *
       * @{
       */

      std::string getName() const
      {
        return name_;
      }
		
      void setName(const std::string& name)
      {
        name_=name;
      }

      /** @} */

      std::vector<std::string> getAllLeavesNames() const;

      void resetNodesId();

      void setBranchLengths(double l);

      /**
       * @brief Multiply all branch lengths by a given factor.
       *
       * @param factor The factor to multiply all branch lengths with.
       */


      void scaleTree(double factor);

      /**
       * @brief Multiply all branch lengths under a Node by a given factor.
       *
       * @param node The node defining the subtree.
       * @param factor The factor to multiply all branch lengths with.
       */

      void scaleTree(std::shared_ptr<PhyloNode> node, double factor);

  };
    
}

#endif
