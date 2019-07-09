#pragma once

#include <string>
#include <vector>
#include <IO/FamiliesFileParser.hpp>
#include <util/enums.hpp>

class DTLRatesVector;


class Routines {
public:
  /**
   * Optimize the DTL rates for the families families. 
   * The result is stored into rates
   */
  static void optimizeRates(bool userDTLRates, 
    const std::string &speciesTreeFile,
    RecModel recModel,
    std::vector<FamiliesFileParser::FamilyInfo> &families,
    bool perSpeciesRates, 
    DTLRatesVector &rates,
    long &sumElapsed);

  /**
   * Infer the reconciliation between the families gene trees and the 
   * species tree, and output them in a file. 
   */
  static void inferReconciliation(
    const std::string &speciesTreeFile,
    std::vector<FamiliesFileParser::FamilyInfo> &families,
    RecModel model,
    DTLRatesVector &rates,
    const std::string &outputDir
    );

  /**
   * Create random trees for families that need one, write them to a file,
   * update the families current gene tree, and return true if there was at
   * least one random tree to generate
   */
  static bool createRandomTrees(const std::string &geneRaxOutputDir, 
      std::vector<FamiliesFileParser::FamilyInfo> &families);
  

  /**
   * Read the family stats files and sum the sequence and reconciliation likelihoods
   */
  static void gatherLikelihoods(std::vector<FamiliesFileParser::FamilyInfo> &families,
    double &totalLibpllLL,
    double &totalRecLL);


};
