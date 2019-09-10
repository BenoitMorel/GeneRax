#pragma once

#include <string>
#include <vector>
#include <IO/FamiliesFileParser.hpp>
#include <util/enums.hpp>

class Parameters;


class Routines {
public:
  /**
   * Optimize the DTL rates for the families families. 
   * The result is stored into rates
   */
  static void optimizeRates(bool userDTLRates, 
    const std::string &speciesTreeFile,
    RecModel recModel,
    Families &families,
    bool perSpeciesRates, 
    Parameters &rates,
    long &sumElapsed);

  /**
   * Infer the reconciliation between the families gene trees and the 
   * species tree, and output them in a file. 
   */
  static void inferReconciliation(
    const std::string &speciesTreeFile,
    Families &families,
    RecModel model,
    Parameters &rates,
    const std::string &outputDir
    );

  /**
   * Create random trees for families that need one, write them to a file,
   * update the families current gene tree, and return true if there was at
   * least one random tree to generate
   */
  static bool createRandomTrees(const std::string &geneRaxOutputDir, 
      Families &families);
  

  /**
   * Read the family stats files and sum the sequence and reconciliation likelihoods
   */
  static void gatherLikelihoods(Families &families,
    double &totalLibpllLL,
    double &totalRecLL);


};
