
#include "Routines.hpp"
#include <sstream>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <trees/SpeciesTree.hpp>
#include <parallelization/ParallelContext.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <optimizers/SpeciesTreeOptimizer.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <IO/FileSystem.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <trees/PLLRootedTree.hpp>
#include <maths/ModelParameters.hpp>
#include <routines/scheduled_routines/RaxmlMaster.hpp>
#include <routines/scheduled_routines/GeneRaxMaster.hpp>
#include <maths/Random.hpp>
#include <trees/PLLRootedTree.hpp>
#include <NJ/MiniNJ.hpp>
#include <NJ/Cherry.hpp>
#include <NJ/CherryPro.hpp>
  
std::unique_ptr<PLLRootedTree> 
Routines::computeInitialSpeciesTree(Families &families,
    const std::string globalOutputDir,
    SpeciesTreeAlgorithm algo)
{
  std::string cladeOutput = FileSystem::joinPaths(
      globalOutputDir, "cladesSpeciesTree");
  switch (algo) {
  case SpeciesTreeAlgorithm::MiniNJ:
    return MiniNJ::runMiniNJ(families); 
  case SpeciesTreeAlgorithm::NJst:
    return MiniNJ::runNJst(families); 
  case SpeciesTreeAlgorithm::WMinNJ:
    return MiniNJ::runWMinNJ(families); 
  case SpeciesTreeAlgorithm::Ustar:
    return MiniNJ::runUstar(families); 
  case SpeciesTreeAlgorithm::Cherry:
    return Cherry::geneTreeCherry(families); 
  case SpeciesTreeAlgorithm::CherryPro:
    return CherryPro::geneTreeCherryPro(families); 
  case SpeciesTreeAlgorithm::Random:
    return std::make_unique<PLLRootedTree>(
        std::make_unique<SpeciesTree>(families)->getTree().getNewickString(), 
        false);
  case SpeciesTreeAlgorithm::User:
    assert(false);
  }
  return nullptr;
}

void Routines::runRaxmlOptimization(Families &families,
    const std::string &output,
    const std::string &execPath,
    unsigned int iteration,
    bool splitImplem,
    long &sumElapsedSec)
{
  RaxmlMaster::runRaxmlOptimization(families, output,
      execPath,
      iteration,
      splitImplem, 
      sumElapsedSec);
}
  
void Routines::optimizeGeneTrees(Families &families,
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
    bool inPlace)
{
  GeneRaxMaster::optimizeGeneTrees(families,
      recModelInfo,
      rates,
      output,
      resultName,
      execPath,
      speciesTreePath,
      reconciliationOpt,
      madRooting,
      supportThreshold,
      recWeight,
      enableRec,
      enableLibpll,
      sprRadius,
      iteration,
      schedulerSplitImplem,
      elapsed,
      inPlace);
}

void Routines::exportPerSpeciesRates(const std::string &speciesTreeFile,
    Parameters &rates,
    const RecModelInfo &recModelInfo,
    const std::string &outputFile)
{
  Logger::info << "Exporting per-species rates into " << outputFile << std::endl;
  ParallelOfstream os(outputFile);
  PLLRootedTree speciesTree(speciesTreeFile);
  auto freeParameters = recModelInfo.modelFreeParameters();
  auto speciesNodesNumber = speciesTree.getNodesNumber();
  assert(speciesNodesNumber * freeParameters == rates.dimensions());
  for (auto node: speciesTree.getNodes()) {
    auto e = node->node_index;
    os << node->label;
    for (unsigned int i = 0; i < freeParameters; ++i) {
      os << " " << rates[e * freeParameters + i];
    }
    os << "\n";
  }
}

void Routines::optimizeRates(bool userDTLRates, 
    const std::string &speciesTreeFile,
    const RecModelInfo &recModelInfo,
    Families &families,
    bool perSpeciesRates, 
    Parameters &rates,
    long &sumElapsed) 
{
  if (userDTLRates) {
    return;
  }
  auto start = Logger::getElapsedSec();
  PerCoreGeneTrees geneTrees(families);
  PLLRootedTree speciesTree(speciesTreeFile);
  speciesTree.ensureUniqueLabels();
  PerCoreEvaluations evaluations;
  buildEvaluations(geneTrees, speciesTree, recModelInfo, evaluations);
  if (perSpeciesRates) {
    assert(speciesTree.areNodeIndicesParallelConsistent());
    rates = DTLOptimizer::optimizeParametersPerSpecies(evaluations, speciesTree.getNodesNumber());
  } else {
    rates = DTLOptimizer::optimizeParametersGlobalDTL(evaluations);
  }
  ParallelContext::barrier(); 
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
}


static std::string getSpeciesEventCountFile(const std::string &outputDir, const std::string &familyName)
{
  return FileSystem::joinPaths(outputDir, 
      FileSystem::joinPaths("reconciliations", familyName + "_speciesEventCounts.txt"));
}


static std::string getTransfersFile(const std::string &outputDir, const std::string &familyName, int sample = -1)
{
  auto res = FileSystem::joinPaths(outputDir, FileSystem::joinPaths("reconciliations", familyName));
  if (sample >= 0) {
    res += std::string("_") + std::to_string(sample);
  }
  res += std::string("_transfers.txt");
  return res;
}
void Routines::inferAndGetReconciliationScenarios(
    PLLRootedTree &speciesTree,
    const PerCoreGeneTrees &geneTrees,
    const ModelParameters &initialModelRates,
    unsigned int reconciliationSamples,
    bool optimizeRates,
    std::vector<Scenario> &scenarios)
{
  // initialization
  auto consistentSeed = Random::getInt();
  auto modelParameters = initialModelRates;
  ParallelContext::barrier();
  // pre-fill scenarios array
  const unsigned int scenariosNumber = geneTrees.getTrees().size() * 
    (std::max(1u, reconciliationSamples));
  scenarios = std::vector<Scenario>(scenariosNumber);
  // initialize evaluations, and optimize if requested
  PerCoreEvaluations evaluations;
  evaluations.resize(geneTrees.getTrees().size());
  for (unsigned int i = 0; i  < geneTrees.getTrees().size(); ++i) {
    auto &tree = geneTrees.getTrees()[i];
    evaluations[i] = std::make_shared<ReconciliationEvaluation>(speciesTree, *tree.geneTree, tree.mapping, modelParameters.info);
    evaluations[i]->setRates(modelParameters.getRates(i));
  }
  if (optimizeRates) {
    Logger::timed << "Optimizing DTL rates before reconciliation..." 
      << std::endl;
    bool optimizeFromStartingParameters = true;
    modelParameters = DTLOptimizer::optimizeModelParameters(
        evaluations, 
        optimizeFromStartingParameters, 
        modelParameters);
    Logger::timed << "done" << std::endl;
  }
  ParallelContext::barrier();
  // infer the scenarios!
  for (unsigned int i = 0; i  < geneTrees.getTrees().size(); ++i) {
    if (reconciliationSamples < 1) {
      evaluations[i]->inferMLScenario(scenarios[i]);
    } else {
      for (unsigned int sample = 0; sample < reconciliationSamples; ++sample) {
        auto index = sample * geneTrees.getTrees().size() + i;
        evaluations[i]->inferMLScenario(scenarios[index], true);
      }
    }
  }
  // restore the seed to a consistent state
  Random::setSeed(consistentSeed);
  ParallelContext::barrier();
}

  

void Routines::inferReconciliation(
    const std::string &speciesTreeFile,
    Families &families,
    const ModelParameters &initialModelRates,
    const std::string &outputDir,
    bool bestReconciliation,
    unsigned int reconciliationSamples,
    bool optimizeRates
    )
{
  PerCoreGeneTrees geneTrees(families);
  std::string reconciliationsDir = FileSystem::joinPaths(outputDir, "reconciliations");
  FileSystem::mkdir(reconciliationsDir, true);
  PLLRootedTree speciesTree(speciesTreeFile);
  if (bestReconciliation) {
    std::vector<Scenario> scenarios;
    inferAndGetReconciliationScenarios(speciesTree, 
        geneTrees, 
        initialModelRates,
        0,
        optimizeRates, 
        scenarios);
    assert(scenarios.size() == geneTrees.getTrees().size());
    for (unsigned int i = 0; i  < geneTrees.getTrees().size(); ++i) {
      auto &tree = geneTrees.getTrees()[i];
      std::string eventCountsFile = FileSystem::joinPaths(reconciliationsDir, tree.name + "_eventCounts.txt");
      std::string speciesEventCountsFile = getSpeciesEventCountFile(outputDir, tree.name);
      std::string transfersFile = getTransfersFile(outputDir, tree.name);
      std::string orthoGroupFile = FileSystem::joinPaths(reconciliationsDir, tree.name + "_orthogroups.txt");
      std::string allOrthoGroupFile = FileSystem::joinPaths(reconciliationsDir, tree.name + "_orthogroups_all.txt");
      std::string treeWithEventsFileNHX = FileSystem::joinPaths(reconciliationsDir, tree.name + "_reconciliated.nhx");
      std::string treeWithEventsFileRecPhyloXML = FileSystem::joinPaths(reconciliationsDir, 
          tree.name + "_reconciliated.xml");
      std::string treeWithEventsFileNewickEvents = FileSystem::joinPaths(
          reconciliationsDir, 
          tree.name + "_events.newick");
      auto &scenario = scenarios[i];
      scenario.saveEventsCounts(eventCountsFile, false);
      scenario.savePerSpeciesEventsCounts(speciesEventCountsFile, false);
      scenario.saveReconciliation(treeWithEventsFileRecPhyloXML, ReconciliationFormat::RecPhyloXML, false);
      scenario.saveReconciliation(treeWithEventsFileNHX, ReconciliationFormat::NHX, false);
      scenario.saveLargestOrthoGroup(orthoGroupFile, false);
      scenario.saveAllOrthoGroups(allOrthoGroupFile, false);
      scenario.saveReconciliation(treeWithEventsFileNewickEvents, 
          ReconciliationFormat::NewickEvents, false);
      scenario.saveReconciliation(treeWithEventsFileNHX, ReconciliationFormat::NHX, false);
      scenario.saveTransfers(transfersFile, false);
    }
  }
  if (reconciliationSamples) {
    std::vector<Scenario> scenarios;
    inferAndGetReconciliationScenarios(speciesTree, 
        geneTrees, 
        initialModelRates, 
        reconciliationSamples, 
        optimizeRates, 
        scenarios);
    assert(scenarios.size() == geneTrees.getTrees().size() * reconciliationSamples);
    for (unsigned int i = 0; i  < geneTrees.getTrees().size(); ++i) {
      auto &tree = geneTrees.getTrees()[i];
      std::string nhxSamples = FileSystem::joinPaths(reconciliationsDir, tree.name + "_samples.nhx");
      ParallelOfstream nhxOs(nhxSamples, false);
      for (unsigned int sample = 0; sample < reconciliationSamples; ++sample) {
        auto scenarioIndex = sample * geneTrees.getTrees().size() + i;
        auto &scenario = scenarios[scenarioIndex];
        std::string transfersFile = getTransfersFile(outputDir, tree.name, sample);
        scenario.saveReconciliation(nhxOs, ReconciliationFormat::NHX);
        scenario.saveTransfers(transfersFile, false);
        scenario.resetBlackList();
        nhxOs << "\n";
      }
    }
  }
  ParallelContext::barrier();
  {
    bool forceTransfers = false;
    PerSpeciesEvents events;
    getPerSpeciesEvents(speciesTree,
      geneTrees,
      initialModelRates,
      reconciliationSamples,
      events,
      forceTransfers);
    PLLRootedTree speciesTree(speciesTreeFile);
    ParallelOfstream os(FileSystem::joinPaths(outputDir, "per_species_event_counts.txt"));
    os << "#S #SL #D #T #TL" << std::endl;
    for (unsigned int e = 0; e < events.events.size(); ++e) {
      std::string label(speciesTree.getNode(e)->label);
      auto &eventCount = events.events[e];
      os << label << " ";
      os << "S=" << eventCount.SCount + eventCount.LeafCount << " ";
      os << "SL=" << eventCount.SLCount << " ";
      os << "D=" << eventCount.DCount << " ";
      os << "T=" << eventCount.TCount << " ";
      os << "TL=" << eventCount.TLCount << std::endl;
    }
  }
}
  
void Routines::computeSuperMatrixFromOrthoGroups(
      const std::string &speciesTreeFile,
      Families &families,
      const std::string &outputDir,
      const std::string &outputFasta,
      bool largestOnly,
      bool masterOnly)
{
  auto savedSeed = Random::getInt(); // for some reason, parsing
                          // the Model calls rand, so we
                          // have so ensure seed consistency
  Random::setSeed(savedSeed);
  if (masterOnly && ParallelContext::getRank() != 0) {
    return;
  }
  PLLRootedTree speciesTree(speciesTreeFile);
  const std::unordered_set<std::string> speciesSet(
      speciesTree.getLabels(true));
  std::string reconciliationsDir = FileSystem::joinPaths(outputDir, "reconciliations");
  std::unordered_map<std::string, std::string> superMatrix;
  for (auto &species: speciesSet) {
    superMatrix.insert({species, std::string()});
  }
  unsigned int offset = 0;
  unsigned int currentSize = 0;
  std::ofstream partitionOs(outputFasta + ".part");
  for (auto &family: families) {
    std::string orthoGroupFile = FileSystem::joinPaths(reconciliationsDir, family.name);
    if (largestOnly) {
      orthoGroupFile +=  "_orthogroups.txt";
    } else {
      orthoGroupFile +=  "_orthogroups_all.txt";
    }
    OrthoGroups orthoGroups;
    parseOrthoGroups(orthoGroupFile, orthoGroups);
    for (auto &orthoGroup: orthoGroups) {
      if (orthoGroup->size() < 4) {
        continue;
      }
      auto model = LibpllParsers::getModel(family.libpllModel);
      PLLSequencePtrs sequences;
      unsigned int *weights = nullptr;
      LibpllParsers::parseMSA(family.alignmentFile, 
        model->charmap(),
        sequences,
        weights);
      GeneSpeciesMapping mapping;
      mapping.fill(family.mappingFile, family.startingGeneTree);
      for (auto &sequence: sequences) {
        std::string geneLabel(sequence->label);
        if (orthoGroup->find(geneLabel) != orthoGroup->end()) {
          // add the sequence to the supermatrix
          superMatrix[mapping.getSpecies(geneLabel)] += std::string(sequence->seq);
          offset = superMatrix[mapping.getSpecies(geneLabel)].size();
          currentSize = sequence->len;
        }
      }
      partitionOs << model->name() << ", " << family.name;
      partitionOs << " = " << offset - currentSize + 1 << "-" << offset << std::endl;
      std::string gaps(currentSize, '-');
      for (auto &superPair: superMatrix) {
        auto &superSequence = superPair.second;
        if (superSequence.size() != offset) {
          superSequence += gaps;
          assert(superSequence.size() == offset);
        }
      }
      free(weights);
    }
  }
  LibpllParsers::writeSuperMatrixFasta(superMatrix, outputFasta);
  Random::setSeed(savedSeed);
}


bool Routines::createRandomTrees(const std::string &geneRaxOutputDir, Families &families)
{
  std::string startingTreesDir = FileSystem::joinPaths(geneRaxOutputDir, "startingTrees");
  bool startingTreesDirCreated = false;
  auto consistentSeed = Random::getInt();
  for (auto &family: families) {
    if (family.startingGeneTree == "__random__") {
        if (!startingTreesDirCreated) {
          FileSystem::mkdir(startingTreesDir, true);
          startingTreesDirCreated = true;
        } 
        family.startingGeneTree = FileSystem::joinPaths(geneRaxOutputDir, "startingTrees");
        family.startingGeneTree = FileSystem::joinPaths(family.startingGeneTree, family.name + ".newick");
        if (ParallelContext::getRank() == 0) {
          LibpllEvaluation::createAndSaveRandomTree(family.alignmentFile, family.libpllModel, family.startingGeneTree);
        }
    }
  }
  Random::setSeed(consistentSeed);
  ParallelContext::barrier();
  return startingTreesDirCreated;
}

void Routines::gatherLikelihoods(Families &families,
    double &totalLibpllLL,
    double &totalRecLL)
{
  ParallelContext::barrier();
  totalRecLL = 0.0;
  totalLibpllLL = 0.0;
  unsigned int familiesNumber = static_cast<unsigned int>(families.size());
  for (auto i = ParallelContext::getBegin(familiesNumber); i < ParallelContext::getEnd(familiesNumber); ++i) {
    auto &family = families[i];
    std::ifstream is(family.statsFile);
    double libpllLL = 0.0;
    double recLL = 0.0;
    is >> libpllLL;
    is >> recLL;
    totalRecLL += recLL;
    totalLibpllLL += libpllLL;
  }
  ParallelContext::sumDouble(totalRecLL);
  ParallelContext::sumDouble(totalLibpllLL);
}
  

static const std::string keyDelimiter("-_-");


void Routines::getLabelsFromTransferKey(const std::string &key, std::string &label1, std::string &label2)
{
  auto pos = key.find(keyDelimiter);
  label1 = key.substr(0, pos);
  label2 = key.substr(pos + keyDelimiter.size());
}



void Routines::getPerSpeciesEvents(PLLRootedTree &speciesTree,
  PerCoreGeneTrees &geneTrees,
  const ModelParameters &modelParameters,
  unsigned int reconciliationSamples,
  PerSpeciesEvents &events,
  bool forceTransfers)
{
  events = PerSpeciesEvents(speciesTree.getNodesNumber());
  std::vector<Scenario> scenarios;
  bool optimizeRates = false;
  ModelParameters transfersModelParameter;
  if (Enums::accountsForTransfers(modelParameters.info.model)
      || !forceTransfers) {
    transfersModelParameter = modelParameters;
  } else {
    // the current model does not account for transfers
    // to infer transfers, we use a model that does
    Parameters transfersParameters(0.2, 0.2, 0.2);
    RecModelInfo recModelInfo = modelParameters.info;
    recModelInfo.model = RecModel::UndatedDTL;
    recModelInfo.perFamilyRates = false;
    transfersModelParameter = ModelParameters(transfersParameters,
        1,
        recModelInfo);
  }
  transfersModelParameter.info.rootedGeneTree = false;
  inferAndGetReconciliationScenarios(
    speciesTree,
    geneTrees,
    transfersModelParameter,
    reconciliationSamples,
    optimizeRates,
    scenarios);
  for (auto &scenario: scenarios) {
    scenario.gatherReconciliationStatistics(events);
  }
  ParallelContext::barrier();
  events.parallelSum();
}
  

void Routines::getTransfersFrequencies(PLLRootedTree &speciesTree,
    PerCoreGeneTrees &geneTrees,
    const ModelParameters &modelParameters,
    unsigned int reconciliationSamples,
    TransferFrequencies &transferFrequencies)
{
  const bool optimizeRates = false;
  ModelParameters transfersModelParameter;
  if (Enums::accountsForTransfers(modelParameters.info.model)) {
    transfersModelParameter = modelParameters;
  } else {
    // the current model does not account for transfers
    // to infer transfers, we use a model that does
    Parameters transfersParameters(0.2, 0.2, 0.2);
    RecModelInfo recModelInfo = modelParameters.info;
    recModelInfo.model = RecModel::UndatedDTL;
    recModelInfo.perFamilyRates = false;
    transfersModelParameter = ModelParameters(transfersParameters,
        1,
        recModelInfo);
  }
  transfersModelParameter.info.rootedGeneTree = false;
  std::vector<Scenario> scenarios;
  inferAndGetReconciliationScenarios(speciesTree, 
      geneTrees, 
      transfersModelParameter,
      reconciliationSamples,
      optimizeRates, 
      scenarios);
  
  const auto labelToId = speciesTree.getDeterministicLabelToId();
  const auto idToLabel = speciesTree.getDeterministicIdToLabel();
  const unsigned int labelsNumber = idToLabel.size();
  const VectorUint zeros(labelsNumber, 0);
  transferFrequencies.count = MatrixUint(labelsNumber, zeros);
  transferFrequencies.idToLabel = idToLabel;
  unsigned int maxSample = std::max<unsigned int>(1u, reconciliationSamples);
  for (unsigned int sample = 0; sample < maxSample; ++sample) {
    for (unsigned int i = 0; i < geneTrees.getTrees().size(); ++i) {
      auto index = sample * geneTrees.getTrees().size() + i;
      auto &scenario = scenarios[index];
      scenario.countTransfers(labelToId, 
          transferFrequencies.count);
    }
  }
  for (unsigned int i = 0; i < labelsNumber; ++i) {
    ParallelContext::sumVectorUInt(transferFrequencies.count[i]); 
  }
  ParallelContext::barrier();
  assert(ParallelContext::isRandConsistent());
}


void Routines::buildEvaluations(PerCoreGeneTrees &geneTrees, 
    PLLRootedTree &speciesTree, 
    const RecModelInfo &recModelInfo,
    Evaluations &evaluations)
{
  auto &trees = geneTrees.getTrees();
  evaluations.resize(trees.size());
  for (unsigned int i = 0; i < trees.size(); ++i) {
    auto &tree = trees[i];
    evaluations[i] = std::make_shared<ReconciliationEvaluation>(speciesTree, 
        *tree.geneTree, 
        tree.mapping, 
        recModelInfo);
  }
}


void Routines::parseOrthoGroups(const std::string &familyName,
      OrthoGroups &orthoGroups)
{
  std::ifstream is(familyName);
  std::string tmp;
  OrthoGroupPtr currentOrthoGroup = std::make_shared<OrthoGroup>();
  while (is >> tmp) {
    if (tmp == std::string("-")) {
      orthoGroups.push_back(currentOrthoGroup);
      currentOrthoGroup = std::make_shared<OrthoGroup>();
    } else {
      currentOrthoGroup->insert(tmp);
    }
  }
  //orthoGroups.push_back(currentOrthoGroup);
}

