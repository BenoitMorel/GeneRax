//
// File: NexusIOTree.h
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

#ifndef _NEXUSIOTREE_H_
#define _NEXUSIOTREE_H_

#include "IoTree.h"
#include "../Tree/PhyloTree.h"

namespace bpp
{

  /**
   * @brief a simple parser for reading trees from a Nexus file.
   *
   * This reader is not supposed to be a full parser of the Nexus files,
   * but only extract the tree data. Only a basic subset of the options
   * are and will be supported.
   *
   * This format is described in the following paper:
   * Maddison D, Swofford D, and Maddison W (1997), _Syst Biol_ 46(4):590-621
   *
   * @author Julien Dutheil
   */
  class NexusIOTree:
    public virtual AbstractITree,
    public virtual AbstractOTree,
    public virtual AbstractIMultiTree,
    public virtual AbstractOMultiTree
  {
  public:
		
    /**
     * @brief Build a new Nexus tree parser.
     */
    NexusIOTree() {}

    virtual ~NexusIOTree() {}
  
  public:

    /**
     * @name The IOTree interface
     *
     * @{
     */
    const std::string getFormatName() const;
    const std::string getFormatDescription() const;
    /* @} */

    /**
     * @name The ITree interface
     *
     * @{
     */

    PhyloTree* readP(const std::string& path) const
    {
      return AbstractITree::readP(path);
    }
    
    PhyloTree* readP(std::istream& in) const;
    
/** @} */

    /**
     * @name The OTree interface
     *
     * @{
     */
    void write(const PhyloTree& tree, const std::string& path, bool overwrite = true) const
    {
      AbstractOTree::write(tree, path, overwrite);
    }
    void write(const PhyloTree& tree, std::ostream& out) const
    {
      write_(tree, out);
    }
    /** @} */

    /**
     * @name The IMultiTree interface
     *
     * @{
     */

    void read(const std::string& path, std::vector<PhyloTree*>& trees) const
    {
      AbstractIMultiTree::read(path, trees);
    }
    
    void read(std::istream& in, std::vector<PhyloTree*>& trees) const;
/**@}*/

    /**
     * @name The OMultiTree interface
     *
     * @{
     */
    void write(const std::vector<const PhyloTree*>& trees, const std::string& path, bool overwrite = true) const
    {
      AbstractOMultiTree::write(trees, path, overwrite);
    }
    void write(const std::vector<const PhyloTree*>& trees, std::ostream& out) const
    {
      write_(trees, out);
    }
    /** @} */

  protected:

    void write_(const PhyloTree& tree, std::ostream& out) const;

    void write_(const std::vector<const PhyloTree*>& trees, std::ostream& out) const;

  };

} //end of namespace bpp.

#endif  //_NEXUSIOTREE_H_

