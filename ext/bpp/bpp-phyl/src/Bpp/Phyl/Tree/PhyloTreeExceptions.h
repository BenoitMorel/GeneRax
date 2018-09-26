//
// File: PhyloTreeExceptions.h
// Created by: Julien Dutheil
// Created on: Mon Nov  3 17:04:46 2003
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

#ifndef _PHYLOTREEEXCEPTIONS_H_
#define _PHYLOTREEEXCEPTIONS_H_

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>

// From the STL:
#include <string>
#include <memory>

namespace bpp
{
  class PhyloNode;
  class PhyloBranch;
  class PhyloTree;

/**
 * @brief General exception thrown when something is wrong with a particular node.
 */
  class PhyloNodeException :
    public Exception
  {
  protected:
    unsigned int nodeId_;

  public:
    /**
     * @brief Build a new PhyloNodeException.
     * @param text A message to be passed to the exception hierarchy.
     * @param nodeId The id of the node that threw the exception.
     */
    PhyloNodeException(const std::string& text, int nodeId) :
      Exception("PhyloNodeException: " + text + "(id:" + TextTools::toString(nodeId) + ")"),
      nodeId_(nodeId) {}

    virtual ~PhyloNodeException() throw () {}

  public:
    /**
     * @brief Get the id of node that threw the exception.
     *
     * @return The id of the faulty node.
     */
    virtual int getNodeId() const { return nodeId_; }
  };


/**
 * @brief General exception thrown when something is wrong with a particular node.
 */
  class PhyloNodePException :
    public PhyloNodeException
  {
  private:
    const PhyloNode* node_;

  public:
    /**
     * @brief Build a new PhyloNodePException.
     * @param text A message to be passed to the exception hierarchy.
     * @param tree the phyloTree owning the branch
     * @param node A const pointer toward the node that threw the exception.
     */
    PhyloNodePException(const std::string& text, const PhyloTree& tree, const std::shared_ptr<PhyloNode> node);

    /**
     * @brief Build a new PhyloNodePException.
     * @param text A message to be passed to the exception hierarchy.
     * @param node A const pointer toward the node that threw the exception.
     */
    PhyloNodePException(const std::string& text, const PhyloNode* node);

    /**
     * @brief Build a new PhyloNodePException.
     * @param text A message to be passed to the exception hierarchy.
     * @param nodeId The id of the node that threw the exception.
     */
    PhyloNodePException(const std::string& text, int nodeId) :
      PhyloNodeException(text, nodeId), node_(0) {}

    PhyloNodePException(const PhyloNodePException& nex) :
      PhyloNodeException(nex),
      node_(nex.node_)
    {}

    PhyloNodePException& operator=(const PhyloNodePException& nex)
    {
      PhyloNodeException::operator=(nex);
      node_ = nex.node_;
      return *this;
    }

    virtual ~PhyloNodePException() throw () {}

  public:
    /**
     * @brief Get the node that threw the exception.
     *
     * @return A pointer toward the faulty node.
     */
    virtual const PhyloNode* getNode() const { return node_; };
    /**
     * @brief Get the id of node that threw the exception.
     *
     * @return The id of the faulty node.
     */
    virtual int getNodeId() const { return nodeId_; }
  };

/**
 * @brief General exception thrown if a property could not be found.
 */
  class PhyloNodePropertyNotFoundException :
    public PhyloNodePException
  {
  private:
    std::string propertyName_;

  public:
    /**
     * @brief Build a new PropertyNotFoundException (Node).
     * @param text A message to be passed to the exception hierarchy.
     * @param propertyName The name of the property.
     * @param tree the phyloTree owning the node
     * @param node A const pointer toward the node that threw the exception.
     */
    PhyloNodePropertyNotFoundException(const std::string& text, const std::string& propertyName, const PhyloTree& tree, const std::shared_ptr<PhyloNode> node) :
      PhyloNodePException("Property not found: " + propertyName + ". " + text, tree, node),
      propertyName_(propertyName) {}

    /**
     * @brief Build a new PropertyNotFoundException (Node).
     * @param text A message to be passed to the exception hierarchy.
     * @param propertyName The name of the property.
     * @param node A const pointer toward the node that threw the exception.
     */
    PhyloNodePropertyNotFoundException(const std::string& text, const std::string& propertyName, const PhyloNode* node) :
      PhyloNodePException("Property not found: " + propertyName + ". " + text, node),
      propertyName_(propertyName) {}

    /**
     * @brief Build a new PropertyNotFoundException.
     *
     * @param text A message to be passed to the exception hierarchy.
     * @param propertyName The name of the property.
     * @param nodeId The id of the node that threw the exception.
     */
    PhyloNodePropertyNotFoundException(const std::string& text, const std::string& propertyName, int nodeId) :
      PhyloNodePException("Property not found: " + propertyName + ". " + text, nodeId),
      propertyName_(propertyName) {}

    virtual ~PhyloNodePropertyNotFoundException() throw () {}

  public:
    /**
     * @brief Get the name of the property that could not be found.
     *
     * @return The name of the missing property.
     */
    virtual const std::string& getPropertyName() const { return propertyName_; }
  };

/**
 * @brief Exception thrown when something is wrong with a particular node.
 */
  class PhyloNodeNotFoundException :
    public Exception
  {
  private:
    std::string id_;

  public:
    /**
     * @brief Build a new PhyloNodeNotFoundException.
     *
     * @param text A message to be passed to the exception hierarchy.
     * @param id   A string describing the node.
     */
    PhyloNodeNotFoundException(const std::string& text, const std::string& id);

    /**
     * @brief Build a new PhyloNodeNotFoundException.
     *
     * @param text A message to be passed to the exception hierarchy.
     * @param id   A node identifier.
     */
    PhyloNodeNotFoundException(const std::string& text, int id);

    virtual ~PhyloNodeNotFoundException() throw () {}

  public:
    /**
     * @brief Get the node id that threw the exception.
     *
     * @return The id of the node.
     */
    virtual std::string getId() const { return id_; }
  };

  /**
   * @brief General exception thrown when something is wrong with a particular branch.
   */
  class PhyloBranchException :
    public Exception
  {
  protected:
    unsigned int branchId_;

  public:
    /**
     * @brief Build a new PhyloBranchPException.
     * @param text A message to be passed to the exception hierarchy.
     * @param branchId The id of the branch that threw the exception.
     */
    PhyloBranchException(const std::string& text, int branchId) :
      Exception("PhyloBranchException: " + text + "(id:" + TextTools::toString(branchId) + ")"),
      branchId_(branchId) {}

    virtual ~PhyloBranchException() throw () {}

  public:
    /**
     * @brief Get the id of branch that threw the exception.
     *
     * @return The id of the faulty branch.
     */
    virtual int getBranchId() const { return branchId_; }
  };


/**
 * @brief General exception thrown when something is wrong with a particular branch.
 */
  class PhyloBranchPException :
    public PhyloBranchException
  {
  private:
    const PhyloBranch* branch_;

  public:
    /**
     * @brief Build a new PhyloBranchPException.
     * @param text A message to be passed to the exception hierarchy.
     * @param tree the phyloTree owning the branch
     * @param branch A const pointer toward the branch that threw the exception.
     */
    PhyloBranchPException(const std::string& text, const PhyloTree& tree, const std::shared_ptr<PhyloBranch> branch);

    /**
     * @brief Build a new PhyloBranchPException.
     * @param text A message to be passed to the exception hierarchy.
     * @param branch A const pointer toward the branch that threw the exception.
     */
    PhyloBranchPException(const std::string& text, const PhyloBranch* branch);

    /**
     * @brief Build a new PhyloBranchPException.
     * @param text A message to be passed to the exception hierarchy.
     * @param branchId The id of the branch that threw the exception.
     */
    PhyloBranchPException(const std::string& text, int branchId) :
      PhyloBranchException(text, branchId), branch_(0) {}

    PhyloBranchPException(const PhyloBranchPException& nex) :
      PhyloBranchException(nex),
      branch_(nex.branch_)
    {}

    PhyloBranchPException& operator=(const PhyloBranchPException& nex)
    {
      PhyloBranchException::operator=(nex);
      branch_ = nex.branch_;
      return *this;
    }

    virtual ~PhyloBranchPException() throw () {}

  public:
    /**
     * @brief Get the branch that threw the exception.
     *
     * @return A pointer toward the faulty branch.
     */
    virtual const PhyloBranch* getBranch() const { return branch_; };
    /**
     * @brief Get the id of branch that threw the exception.
     *
     * @return The id of the faulty branch.
     */
    virtual int getBranchId() const { return branchId_; }
  };

/**
 * @brief General exception thrown if a property could not be found.
 */
  class PhyloBranchPropertyNotFoundException :
    public PhyloBranchPException
  {
  private:
    std::string propertyName_;

  public:
    /**
     * @brief Build a new PropertyNotFoundException (Branch).
     * @param text A message to be passed to the exception hierarchy.
     * @param propertyName The name of the property.
     * @param tree the phyloTree owning the branch
     * @param branch A const pointer toward the branch that threw the exception.
     */

    PhyloBranchPropertyNotFoundException(const std::string& text, const std::string& propertyName, const PhyloTree& tree, const std::shared_ptr<PhyloBranch> branch) :
      PhyloBranchPException("Property not found: " + propertyName + ". " + text, tree, branch),
      propertyName_(propertyName) {}

    /**
     * @brief Build a new PropertyNotFoundException (Branch).
     * @param text A message to be passed to the exception hierarchy.
     * @param propertyName The name of the property.
     * @param branch A const pointer toward the branch that threw the exception.
     */
    
    PhyloBranchPropertyNotFoundException(const std::string& text, const std::string& propertyName, const PhyloBranch* branch) :
      PhyloBranchPException("Property not found: " + propertyName + ". " + text, branch),
      propertyName_(propertyName) {}

    /**
     * @brief Build a new PropertyNotFoundException.
     *
     * @param text A message to be passed to the exception hierarchy.
     * @param propertyName The name of the property.
     * @param branchId The id of the branch that threw the exception.
     */
    PhyloBranchPropertyNotFoundException(const std::string& text, const std::string& propertyName, int branchId) :
      PhyloBranchPException("Property not found: " + propertyName + ". " + text, branchId),
      propertyName_(propertyName) {}

    virtual ~PhyloBranchPropertyNotFoundException() throw () {}

  public:
    /**
     * @brief Get the name of the property that could not be found.
     *
     * @return The name of the missing property.
     */
    virtual const std::string& getPropertyName() const { return propertyName_; }
  };

/**
 * @brief Exception thrown when something is wrong with a particular branch.
 */
  class PhyloBranchNotFoundException :
    public Exception
  {
  private:
    std::string id_;

  public:
    /**
     * @brief Build a new PhyloBranchNotFoundException.
     *
     * @param text A message to be passed to the exception hierarchy.
     * @param id   A string describing the branch.
     */
    PhyloBranchNotFoundException(const std::string& text, const std::string& id);

    /**
     * @brief Build a new PhyloBranchNotFoundException.
     *
     * @param text A message to be passed to the exception hierarchy.
     * @param id   A branch identifier.
     */
    PhyloBranchNotFoundException(const std::string& text, int id);

    virtual ~PhyloBranchNotFoundException() throw () {}

  public:
    /**
     * @brief Get the branch id that threw the exception.
     *
     * @return The id of the branch.
     */
    virtual std::string getId() const { return id_; }
  };

/**
 * @brief General exception thrown when something wrong happened in a tree.
 */
  class PhyloTreeException :
    public Exception
  {
  private:
    const PhyloTree* tree_;

  public:
    /**
     * @brief Build a new PhyloTreeException.
     *
     * @param text A message to be passed to the exception hierarchy.
     * @param tree A const pointer toward the tree that threw the exception.
     */
    PhyloTreeException(const std::string& text, const PhyloTree* tree = 0);

    PhyloTreeException(const PhyloTreeException& tex) :
      Exception(tex),
      tree_(tex.tree_)
    {}

    PhyloTreeException& operator=(const PhyloTreeException& tex)
    {
      Exception::operator=(tex);
      tree_ = tex.tree_;
      return *this;
    }

    virtual ~PhyloTreeException() throw () {}

  public:
    /**
     * @brief Get the tree that threw the exception.
     *
     * @return The faulty tree
     */
    virtual const PhyloTree* getTree() const { return tree_; }
  };

/**
 * @brief Exception thrown when a tree is expected to be rooted.
 */
  class UnrootedPhyloTreeException :
    public PhyloTreeException
  {
  public:
    /**
     * @brief Build a new UnrootedPhyloTreeException.
     *
     * @param text A message to be passed to the exception hierarchy.
     * @param tree A const pointer toward the tree that threw the exception.
     */
    UnrootedPhyloTreeException(const std::string& text, const PhyloTree* tree = 0);

    virtual ~UnrootedPhyloTreeException() throw () {}
  };

} // end of namespace bpp.

#endif  // _PHYLOTREEEXCEPTIONS_H_

