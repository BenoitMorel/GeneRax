#pragma once

#include <string>
#include <vector>
#include <IO/FamiliesFileParser.hpp>


class RaxmlMaster {
public:
  RaxmlMaster() = delete;
  /*
   *  Schedule gene tree inference using
   *  sequences only, with raxml-ng algorithm.
   *  @param families Families descriptions
   *  @param output GeneRax run output directory
   *  @param execPath GeneRax executable
   *  @param iteration unique ID for this call
   *                   will be used to create a directory
   *  @param splitImpl use the MPIScheduler split implementation 
   *                   (or the fork)
   *  @param sumElapsedSec will be incremented by the number of
   *                       seconds spent in this call
   */
  static void runRaxmlOptimization(Families &families,
    const std::string &output,
    const std::string &execPath,
    unsigned int iteration,
    bool splitImplem,
    long &sumElapsedSec);
};
