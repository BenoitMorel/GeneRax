#pragma once

#include <vector>
#include <string>
#include <IO/FamiliesFileParser.hpp>
#include <util/enums.hpp>

class Parameters;

class GeneRaxMaster {
public:
  GeneRaxMaster() = delete;

  static void optimizeGeneTrees(Families &families,
    RecModel recModel,
    Parameters &rates,
    const std::string &output, 
    const std::string &resultName, 
    const std::string &execPath, 
    const std::string &speciesTreePath,
    RecOpt reconciliationOpt,
    bool perFamilyDTLRates,
    bool rootedGeneTree,
    double recWeight,
    bool enableRec,
    bool enableLibpll,
    int sprRadius,
    int iteration,
    bool schedulerSplitImplem,
    long &elapsed,
    bool inPlace = false); 
};
