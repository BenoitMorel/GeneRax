#pragma once

#include <vector>
#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>




class Cherry {
public:
  Cherry() = delete;

  /**
   *
   */
  static std::unique_ptr<PLLRootedTree> geneTreeCherry(const Families &families);
};

