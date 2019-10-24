#pragma once

#include <string>
#include <vector>
#include <IO/FamiliesFileParser.hpp>
#include <util/enums.hpp>
#include <unordered_map>

class Parameters;


typedef std::unordered_map<std::string, unsigned int> TransferFrequencies;



class Routines {

public:
  Routines() = delete;
  
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

  static void optimizeSpeciesRatesEmpirical(const std::string &speciesTreeFile,
    RecModel recModel,
    Families &families,
    Parameters &rates,
    const std::string &outputDir,
    long &sumElapsed);
  
  static void getTransfersFrequencies(const std::string &speciesTreeFile,
    RecModel recModel,
    Families &families,
    const Parameters &rates,
    TransferFrequencies &frequencies,
    const std::string &outputDir);
  
  static void getParametersFromTransferFrequencies(const std::string &speciesTreeFile,
      const TransferFrequencies &frequencies, 
      Parameters &parameters);

  /**
   * Infer the reconciliation between the families gene trees and the 
   * species tree, and output them in a file. 
   */
  static void inferReconciliation(
    const std::string &speciesTreeFile,
    Families &families,
    RecModel model,
    const Parameters &rates,
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
