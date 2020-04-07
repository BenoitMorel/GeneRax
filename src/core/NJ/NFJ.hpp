#pragma once

#include <vector>
#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>




class NFJ {
public:
  NFJ() = delete;

  /**
   *
   */
  static std::unique_ptr<PLLRootedTree> geneTreeNFJ(const Families &families);
};

