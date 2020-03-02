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
    double supportThreshold,
    double recWeight,
    bool enableRec,
    bool enableLibpll,
    unsigned int sprRadius,
    unsigned int iteration,
    bool schedulerSplitImplem,
    long &elapsed,
    bool inPlace = false); 
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
