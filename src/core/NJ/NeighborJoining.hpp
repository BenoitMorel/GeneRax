#pragma once

#include <families/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>

/*
 * Naive NJ implementation. Could be greatly optimized if needed!!!!
 */
class NeighborJoining {
public:
  static std::unique_ptr<PLLRootedTree> countProfileNJ(const Families &families);  
};
