#pragma once

class GeneRaxSlave {
public:
  GeneRaxSlave() = delete;
  static int optimizeGeneTreesMain(int argc, char** argv, void* comm);
};

