#pragma once

#include <families/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>

class NeighborJoining {
public:
  static std::unique_ptr<PLLRootedTree> countProfileNJ(const Families &families);  
};
