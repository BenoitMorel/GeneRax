//
// File: Nhx.h
// Created by: Bastien Boussau
// Created on: Tu Oct 19 10:24:03 2010
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#ifndef _NHX_H_
#define _NHX_H_

#include "IoTree.h"
#include "../Tree/PhyloTree.h"

//From the STL:
#include <set>

namespace bpp
{

/**
 * @brief The so-called 'Nhx - New Hampshire eXtended' parenthetic format. 
 *
 * See http://www.phylosoft.org/Nhx/ for details.
 *
 * Branch lengths and node annotations are supported:
 *
 * ex:
 * <code>
 * ((A:0.1[&&NHX:S=human], B:0.15[&&NHX:S=chimp]):0.2[&&NHX:B=90:D=N:S=primates], C:0.27[&&NHX:S=mouse]);
 * </code>
 *
 * Where e.g. "S=human" means that the sequence A comes from species "human", "B=90" stands for a support value, 
 * and "D=N" means that there was no duplication at the node. Other tags are allowed, see http://www.phylosoft.org/NHX/.
 * By default, comments between "[" and "]" are removed, unless the opening bracket is followed by "&&NHX".
 *
 * Code example:
 * @code
 * #include <Phyl/Nhx.h>
 * #include <Phyl/Tree/Tree.h>
 * 
 * Nhx * NhxReader = new Nhx();
 * try {
 *   Tree * tree = nhxReader->read("MyTestTree.dnd"); // Tree in file MyTestTree.dnd
 *   cout << "Tree has " << tree->getNumberOfLeaves() << " leaves." << endl;
 * } catch (Exception e) {
 *  cout << "Error when reading tree." << endl;
 * }
 * delete tree;
 * delete nhxReader;
 * @endcode
 *
 * All node annotations are stored as node properties, with type bppString for all properties except for support values, where a Number is used.
 *
 */
  class Nhx:
    public AbstractITree,
    public AbstractOTree,
    public AbstractIMultiTree,
    public AbstractOMultiTree
  {
  public:
    struct Property
    {
    public:
      /**
       * @brief The name of the property, which will be used in parsed trees.
       */
      std::string name;
      /**
       * @brief The tag of the property, as it will be found in the tree file.
       */
      std::string tag;
      /**
       * @brief Tells if the property is a branch property instead of a node property.
       */
      bool onBranch;
      /**
       * @brief The type of the property. 0 is string, 1 is integer, 2 is double, 3 is boolean.
       */
      short type;

    public:
      Property(const std::string& pptName, const std::string& pptTag, bool pptOnBranch = false, short pptType = 0):
        name(pptName), tag(pptTag), onBranch(pptOnBranch), type(pptType) {}

      bool operator<(const Property& ppt) const {
        return (name < ppt.name);
      }

    };

  private:
    std::set<Property> supportedProperties_;
    bool useTagsAsPropertyNames_;
    mutable bool hasIds_;

  public:
    
    /**
     * @brief Build a new Nhx reader/writer.
     *
     * Comments between hooks ('[' ']') are ignored.
     *
     * @param useTagsAsPptNames Tells if the NHX tag should be used as a property name in the parsed tree.
     */
    Nhx(bool useTagsAsPptNames = true);
    virtual ~Nhx() {}
  
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

    PhyloTree* parenthesisToPhyloTree(const std::string& description) const;

    std::string treeToParenthesis(const PhyloTree& tree) const;

    void registerProperty(const Property& property) {
      supportedProperties_.insert(property);
    }

    /**
     * @brief Convert property names from tag to names.
     *
     * If a tree has been parsed using useTagsAsPropertyNames=true,
     * this method allows to convert the tree as is it was parsed using
     * the option set to false.
     *
     * @param node The root node of the subtree to convert.
     */
    void changeTagsToNames(PhyloTree& tree, std::shared_ptr<PhyloNode> node) const;

    /**
     * @brief Convert property names from names to tags.
     *
     * If a tree has been parsed using useTagsAsPropertyNames=false,
     * this method allows to convert the tree as is it was parsed using
     * the option set to true.
     *
     * @param node The root node of the subtree to convert.
     */
    void changeNamesToTags(PhyloTree& tree, std::shared_ptr<PhyloNode> node) const;

    void useTagsAsPropertyNames(bool yn) { useTagsAsPropertyNames_ = yn; }
    bool useTagsAsPropertyNames() const { return useTagsAsPropertyNames_; }

  protected:

    void write_(const PhyloTree& tree, std::ostream& out) const;
    
    void write_(const std::vector<const PhyloTree*>& trees, std::ostream& out) const;

    AbstractITree::Element getElement(const std::string& elt) const;

  private:

    std::shared_ptr<PhyloNode> parenthesisToNode(PhyloTree& tree, std::shared_ptr<PhyloNode> father, const std::string& description) const;

public:

    std::string propertiesToParenthesis(const PhyloTree& tree, const std::shared_ptr<PhyloNode> node) const;

protected:

    std::string nodeToParenthesis(const PhyloTree& tree, const std::shared_ptr<PhyloNode> node) const;

    bool setNodeProperties(PhyloTree& tree, std::shared_ptr<PhyloNode> node, const std::string properties) const;

  protected:
    static std::string propertyToString_(const Clonable* pptObject, short type);
    static Clonable* stringToProperty_(const std::string& pptDesc, short type);
  };

} //end of namespace bpp.

#endif  //_NHX_H_

