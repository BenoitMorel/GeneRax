#pragma once

#include <vector>
#include <string>
#include <IO/FamiliesFileParser.hpp>
#include <util/enums.hpp>

class DTLRatesVector;

class GeneTreeSearchMaster {
public:
  static void optimizeGeneTrees(Families &families,
    RecModel recModel,
    DTLRatesVector &rates,
    const std::string &output, 
    const std::string &execPath, 
    const std::string &speciesTreePath,
    RecOpt reconciliationOpt,
    bool rootedGeneTree,
    double recWeight,
    bool enableRec,
    int sprRadius,
    int iteration,
    bool schedulerSplitImplem,
    long &sumElapsed); 
};
