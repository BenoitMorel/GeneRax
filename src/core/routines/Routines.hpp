#pragma once

#include <string>
#include <vector>
#include <IO/FamiliesFileParser.hpp>
#include <util/enums.hpp>
#include <unordered_map>
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <util/Scenario.hpp>


class Parameters;
class ModelParameters;
class PLLRootedTree;
class PerCoreGeneTrees;


class Routines {

public:
  Routines() = delete;


  static std::unique_ptr<PLLRootedTree> computeInitialSpeciesTree(
      Families &family,
      const std::string globalOutputDir,
      SpeciesTreeAlgorithm algo);

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
    const RecModelInfo &recModelInfo,
    Parameters &rates,
    const std::string &output, 
    const std::string &resultName, 
    const std::string &execPath, 
    const std::string &speciesTreePath,
    RecOpt reconciliationOpt,
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
  /**
   * Optimize the DTL rates for the families families. 
   * The result is stored into rates
   */
  static void optimizeRates(bool userDTLRates, 
    const std::string &speciesTreeFile,
    const RecModelInfo &recModelInfo,
    Families &families,
    bool perSpeciesRates, 
    Parameters &rates,
    long &sumElapsed);

  static void getPerSpeciesEvents(PLLRootedTree &speciesTree,
    PerCoreGeneTrees &geneTrees,
    const ModelParameters &modelRates,
    unsigned int reconciliationSamples,
    PerSpeciesEvents &events,
    bool forceTransfers);

  
  static void exportPerSpeciesRates(const std::string &speciesTreeFile,
      Parameters &rates,
      const RecModelInfo &recModelInfo,
      const std::string &outputFile);
  
  static void getTransfersFrequencies(PLLRootedTree &speciesTree,
    PerCoreGeneTrees &geneTrees,
    const ModelParameters &modelRates,
    unsigned int reconciliationSamples,
    TransferFrequencies &frequencies);
  
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
    const ModelParameters &initialModelRates,
    const std::string &outputDir,
    bool bestReconciliation,
    unsigned int reconciliationSamples,
    bool optimizeRates
    );

  /**
   * Infer the reconciliation scenarios between the families gene trees 
   * and the species tree. 
   */
  static void inferAndGetReconciliationScenarios(
    PLLRootedTree &speciesTree,
    const PerCoreGeneTrees &geneTrees,
    const ModelParameters &initialModelRates,
    unsigned int reconciliationSamples, // 0 for ML
    bool optimizeRates,
    std::vector<Scenario> &scenarios);
 
  /**
   *  todobenoit
   */
  static void computeSuperMatrixFromOrthoGroups(
      const std::string &speciesTreeFile,
      Families &families,
      const std::string &outputDir,
      const std::string &outputFasta,
      bool largestOnly,
      bool masterOnly =  true);

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
    const RecModelInfo &recModelInfo,
    Evaluations &evaluations);


private:
  static void parseOrthoGroups(const std::string &familyName,
      OrthoGroups &orthoGroup);

};
