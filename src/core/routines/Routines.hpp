#pragma once

#include <string>
#include <vector>
#include <IO/FamiliesFileParser.hpp>
#include <util/enums.hpp>
#include <unordered_map>
#include <likelihoods/ReconciliationEvaluation.hpp>


class Parameters;
class ModelParameters;
class PLLRootedTree;
class PerCoreGeneTrees;
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
    bool rootedGeneTree,
    bool pruneSpeciesTree,
    Families &families,
    bool perSpeciesRates, 
    Parameters &rates,
    long &sumElapsed);

  static void getTransfersFrequencies(const std::string &speciesTreeFile,
    Families &families,
    const ModelParameters &modelRates,
    TransferFrequencies &frequencies,
    const std::string &outputDir);
  
  static void getLabelsFromTransferKey(const std::string &key, 
      std::string &label1, 
      std::string &label2);
  
  static void getParametersFromTransferFrequencies(const std::string &speciesTreeFile,
      const TransferFrequencies &frequencies, 
      Parameters &parameters);

  /**
   * Infer the reconciliation between the families gene trees and the 
   * species tree, and output them in different files.
   * In addition, perform a stochastich sample of the reconciliations
   */
  static void inferReconciliation(
    const std::string &speciesTreeFile,
    Families &families,
    const ModelParameters &modelRates,
    const std::string &outputDir,
    bool bestReconciliation,
    unsigned int reconciliationSamples,
    bool saveTransfersOnly = false
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

  static void buildEvaluations(PerCoreGeneTrees &geneTrees, 
    PLLRootedTree &speciesTree, 
    RecModel recModel, 
    bool rootedGeneTree, 
    bool pruneSpeciesTree, 
    Evaluations &evaluations);



};
