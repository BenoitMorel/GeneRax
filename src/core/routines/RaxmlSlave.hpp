#pragma once

class RaxmlSlave {
public:
  RaxmlSlave() = delete;
  static int runRaxmlOptimization(int argc, char** argv, void* comm);
};
