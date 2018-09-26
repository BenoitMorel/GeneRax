// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 23/02/17.
//

#ifndef PHYLASOLVER_NODEPROPERTY_H
#define PHYLASOLVER_NODEPROPERTY_H

#include <Bpp/Clonable.h>
#include <memory>
#include <Bpp/Phyl/Tree/PhyloNode.h>
#include "Cost.h"

/*!
 * @class NodeProperty
 * @brief Property of a bpp::PhyloNode.
 * @tparam T
 * @details Used to add properties in PhyloNodes. PhyloNode::getProperty(...) method returns Clonable* elements.
 * To get NodeProperty content, use std::dynamic_cast<...>(...).
 *
 *
 * Example using NodeProperty IsArtificialGene:
 *
 *      auto property = dynamic_cast<IsArtificialGene*>(node->getProperty("isArtificialGene"));
 *
 */
template<typename T>
class NodeProperty: public bpp::Clonable {
private:

protected:
  T value_;

public:
/****************
 * Constructors
 */
  NodeProperty() {}
  NodeProperty(const T& value): value_(value) {}
  NodeProperty(const NodeProperty& model): value_(model.value_) {}


/****************
 * Destructor
 */
  virtual ~NodeProperty() {}


/****************
 * Getters
 */
  T get() const { return value_; };
  virtual NodeProperty * clone() const = 0;
  friend std::ostream& operator<<(std::ostream& os, const NodeProperty& np){
    /// Print NodeProperty content.
    os << np.value_;
    return os;
  }


/****************
 * Setters
 */
  void set(const T& value) { value_ = value; }

};

///////////////////////////////////////////////////////////////

/*!
 * @class IsArtificialGene
 * @brief Node property. Indicates if a gene is artificial (added by the program) or not.
 * @details See NodeProperty.
 */

class IsArtificialGene : public NodeProperty<bool> {
private:

protected:

public:
/****************
 * Constructors
 */
  IsArtificialGene(): NodeProperty() {}
  IsArtificialGene(const bool& value): NodeProperty(value) {}
  IsArtificialGene(const IsArtificialGene& model): NodeProperty(model) {}

/****************
 * Destructor
 */



/****************
 * Getters
 */
  IsArtificialGene * clone() const { return new IsArtificialGene(*this); }


/****************
 * Setters
 */

};

///////////////////////////////////////////////////////////////

/*!
 * @class AssociatedNode
 * @brief Node Property. Indicates the associated node to the current one.
 * @details See NodeProperty.
 */
class AssociatedNode : public NodeProperty<std::shared_ptr<bpp::PhyloNode>> {
private:

protected:

public:
/****************
 * Constructors
 */
  AssociatedNode(): NodeProperty() {}
  AssociatedNode(const std::shared_ptr<bpp::PhyloNode>& value): NodeProperty(value) {}
  AssociatedNode(const AssociatedNode& model): NodeProperty(model) {}

/****************
 * Destructor
 */



/****************
 * Getters
 */
  AssociatedNode* clone() const { return new AssociatedNode(*this); }


/****************
 * Setters
 */

};

///////////////////////////////////////////////////////////////

/*!
 * @class NodeOccurence
 * @brief Node property. Gives the occurence of a node given by PolytomySolver.
 * @details See NodeProperty.
 */

class NodeOccurence : public NodeProperty<std::size_t> {
private:

protected:

public:
/****************
 * Constructors
 */
  NodeOccurence(): NodeProperty(){ this->value_ = 0; }
  NodeOccurence(const std::size_t& value): NodeProperty(value) {}
  NodeOccurence(const NodeOccurence& model): NodeProperty(model) {}

/****************
 * Destructor
 */



/****************
 * Getters
 */
  NodeOccurence * clone() const { return new NodeOccurence(*this); }


/****************
 * Setters
 */

};
#endif //PHYLASOLVER_NODEPROPERTY_H
