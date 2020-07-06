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
    bool madRooting,
    double supportThreshold,
    double recWeight,
    bool enableRec,
    bool enableLibpll,
    unsigned int sprRadius,
    unsigned int iteration,
    bool schedulerSplitImplem,
    long &elapsed,
    bool inPlace = false); 
};
