#pragma once

#include <vector>
#include <IO/Families.hpp>
#include <trees/PLLRootedTree.hpp>
#include <memory>




class CherryPro {
public:
  CherryPro() = delete;

  /**
   *
   */
  static std::unique_ptr<PLLRootedTree> geneTreeCherryPro(const Families &families);
};

