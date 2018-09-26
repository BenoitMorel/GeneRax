//
// File: PhyloNode.h
// Created by: Thomas Bigot
// Created on: Thu Mar 13 12:03:18 2003
//

/*
 *  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)
 * 
 *  This software is a computer program whose purpose is to provide classes
 *  for phylogenetic data analysis.
 * 
 *  This software is governed by the CeCILL  license under French law and
 *  abiding by the rules of distribution of free software.  You can  use, 
 *  modify and/ or redistribute the software under the terms of the CeCILL
 *  license as circulated by CEA, CNRS and INRIA at the following URL
 *  "http://www.cecill.info". 
 * 
 *  As a counterpart to the access to the source code and  rights to copy,
 *  modify and redistribute granted by the license, users are provided only
 *  with a limited warranty  and the software's author,  the holder of the
 *  economic rights,  and the successive licensors  have only  limited
 *  liability. 
 * 
 *  In this respect, the user's attention is drawn to the risks associated
 *  with loading,  using,  modifying and/or developing or reproducing the
 *  software by the user in light of its specific status of free software,
 *  that may mean  that it is complicated to manipulate,  and  that  also
 *  therefore means  that it is reserved for developers  and  experienced
 *  professionals having in-depth computer knowledge. Users are therefore
 *  encouraged to load and test the software's suitability as regards their
 *  requirements in conditions enabling the security of their systems and/or 
 *  data to be ensured and,  more generally, to use and operate it in the 
 *  same conditions as regards security. 
 * 
 *  The fact that you are presently reading this means that you have had
 *  knowledge of the CeCILL license and that you accept its terms.
 */

#ifndef _PHYLONODE_H_
#define _PHYLONODE_H_

#include "PhyloTreeExceptions.h"

#include <Bpp/Clonable.h>
#include <Bpp/Utils/MapTools.h>

#include <memory>

namespace bpp
{
  
  class PhyloNode
  {
  private:
    // a name, if specified
    std::string name_;
    
    // Node properties
    mutable std::map<std::string, Clonable*> properties_;

    
  public:
    /**
     * @brief Build a new void Node object.
     */
    PhyloNode() :
      name_(""),
      properties_()
    {
    }
    
    /**
     * @brief Build a new Node with specified name.
     *
     *
     */
    PhyloNode(const std::string& name):
    name_(name),
    properties_()
    {
    }
    
    /**
     * @brief Copy constructor.
     *
     * @param node The node to copy.
     */
    PhyloNode(const PhyloNode& node);
    
    /**
     * @brief Assignation operator.
     *
     * @param node the node to copy.
     * @return A reference toward this node.
     */

    PhyloNode& operator=(const PhyloNode& node);
    
    PhyloNode* clone() const { return new PhyloNode(*this); }

    /**
     * @brief destructor.
     *
     */

    virtual ~PhyloNode()
    {
      deleteProperties();
    }
    
  public:
    /**
     * @name Name:
     *
     * @{
     */
    
    /**
     * @brief Get the name associated to this node, if there is one,
     * otherwise throw a NodeException.
     *
     * @return The name associated to this node.
     */

    std::string getName() const 
    {
      if (!hasName()) throw PhyloNodePException("Node::getName: no name associated to this node.", this);
      return name_;
    }
    
    /**
     * @brief Give a name or update the name associated to the node.
     *
     * @param name The name to give to the node.
     */

    void setName(const std::string& name)
    {
      name_ = name;
    }
    
    /**
     * @brief Delete the name associated to this node (do nothing if there is no name).
     */

    void deleteName()
    {
      name_ = "";
    }
    
    /**
     * @brief Tell is this node has a name.
     *
     * @return True if name != 0.
     */

    bool hasName() const { return name_ != ""; }
    
    /** @} */
    
    /**
     * @name Node properties:
     *
     * @{
     */
    
    /**
     * @brief Set/add a node property.
     *
     * If no property with the same name is found, the new property will be added to the list.
     * Conversely, the property will be deleted and replaced by the new one.
     * If you want to keep a copy of the old property, consider using the removeProperty function before.
     *
     * @param name The name of the property to set.
     * @param property The property object (will be cloned).
     */

    void setProperty(const std::string& name, const Clonable& property)
    {
      if (hasProperty(name))
        delete properties_[name];
      properties_[name] = property.clone();
    }
    
    Clonable* getProperty(const std::string& name) 
    {
      if (hasProperty(name))
        return properties_[name];
      else
        throw PhyloNodePropertyNotFoundException("", name, this);
    }
    
    const Clonable* getProperty(const std::string& name) const 
    {
      if (hasProperty(name))
        return const_cast<const Clonable*>(properties_[name]);
      else
        throw PhyloNodePropertyNotFoundException("", name, this);
    }
    
    Clonable* removeProperty(const std::string& name) 
    {
      if (hasProperty(name))
      {
        Clonable* removed = properties_[name];
        properties_.erase(name);
        return removed;
      }
      else
        throw PhyloNodePropertyNotFoundException("", name, this);
    }
    
    void deleteProperty(const std::string& name) 
    {
      if (hasProperty(name))
      {
        delete properties_[name];
        properties_.erase(name);
      }
      else
        throw PhyloNodePropertyNotFoundException("", name, this);
    }
    
    /**
     * @brief Remove all node properties.
     *
     * Attached objects will not be deleted.
     */
    void removeProperties()
    {
      properties_.clear();
    }
    
    /**
     * @brief Delete all node properties.
     */
    void deleteProperties()
    {
      for (std::map<std::string, Clonable*>::iterator i = properties_.begin(); i != properties_.end(); i++)
      {
        delete i->second;
      }
      properties_.clear();
    }
    
    bool hasProperty(const std::string& name) const { return properties_.find(name) != properties_.end(); }
    
    std::vector<std::string> getPropertyNames() const { return MapTools::getKeys(properties_); }
    /** @} */
    
  }; //end of class node 
  
  
} //end of namespace bpp.

#else
namespace bpp{ class PhyloNode; }

#endif  //_PHYLONODE_H_
