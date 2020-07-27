
#include "Routines.hpp"
#include <sstream>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <trees/SpeciesTree.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <parallelization/ParallelContext.hpp>
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
      SpeciesTreeAlgorithm algo)
{
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
  bool ok = geneTrees.checkMappings(speciesTreeFile); 
  if (!ok) {
    Logger::info << "INVALID MAPPINGS" << std::endl;
    ParallelContext::abort(42);
  }
  PLLRootedTree speciesTree(speciesTreeFile);
  PerCoreEvaluations evaluations;
  buildEvaluations(geneTrees, speciesTree, recModelInfo, evaluations);
  if (perSpeciesRates) {
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
        auto index = sample + i * reconciliationSamples;
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
    bool optimizeRates,
    bool saveTransfersOnly
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
        reconciliationSamples,
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
      auto &scenario = scenarios[i];
      if (!saveTransfersOnly) {
        scenario.saveEventsCounts(eventCountsFile, false);
        scenario.savePerSpeciesEventsCounts(speciesEventCountsFile, false);
        scenario.saveReconciliation(treeWithEventsFileRecPhyloXML, ReconciliationFormat::RecPhyloXML, false);
        scenario.saveReconciliation(treeWithEventsFileNHX, ReconciliationFormat::NHX, false);
        scenario.saveLargestOrthoGroup(orthoGroupFile, false);
        scenario.saveAllOrthoGroups(allOrthoGroupFile, false);
      }
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
        auto scenarioIndex = sample + i * reconciliationSamples;
        auto &scenario = scenarios[scenarioIndex];
        std::string transfersFile = getTransfersFile(outputDir, tree.name, sample);
        if (!saveTransfersOnly) {
          scenario.saveReconciliation(nhxOs, ReconciliationFormat::NHX);
        }
        scenario.saveTransfers(transfersFile, false);
        scenario.resetBlackList();
        nhxOs << "\n";
      }
    }
  }
  ParallelContext::barrier();
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

static std::string getTransferKey(const std::string &label1, const std::string &label2)
{
  return label1 + keyDelimiter + label2; 
}

void Routines::getLabelsFromTransferKey(const std::string &key, std::string &label1, std::string &label2)
{
  auto pos = key.find(keyDelimiter);
  label1 = key.substr(0, pos);
  label2 = key.substr(pos + keyDelimiter.size());
}


static std::string getLocalTempFile(const std::string &outputDir,
    unsigned int rank)
{
  std::string f = "temp_rank" + std::to_string(rank) + ".txt";
  return FileSystem::joinPaths(outputDir, f); 
}

void mpiMergeTransferFrequencies(TransferFrequencies &frequencies,
    const std::string &outputDir)
{
  std::string tempPath = getLocalTempFile(outputDir, ParallelContext::getRank());
  std::ofstream os(tempPath);
  for (auto &pair: frequencies) {
    os << pair.first << " " << pair.second << std::endl;
  }
  os.close();
  frequencies.clear();
  ParallelContext::barrier();
  for (unsigned int i = 0; i < ParallelContext::getSize(); ++i) {
    std::ifstream is(getLocalTempFile(outputDir, i));
    assert(is.good());
    std::string line;
    while (std::getline(is, line)) {
      std::istringstream iss(line);
      std::string key;
      unsigned int freq;
      iss >> key >> freq;
      if (frequencies.end() == frequencies.find(key)) {
        frequencies.insert({key, freq});
      } else {
        frequencies[key] += freq;
      }
    }
  }
  ParallelContext::barrier();
  remove(tempPath.c_str());
}

void Routines::getPerSpeciesEvents(const std::string &speciesTreeFile,
  Families &families,
  const ModelParameters &modelRates,
  unsigned int reconciliationSamples,
  PerSpeciesEvents &events)
{
  PLLRootedTree speciesTree(speciesTreeFile);
  events = PerSpeciesEvents(speciesTree.getNodesNumber());
  std::vector<Scenario> scenarios;
  bool optimizeRates = false;
  PerCoreGeneTrees geneTrees(families);
  inferAndGetReconciliationScenarios(
    speciesTree,
    geneTrees,
    modelRates,
    reconciliationSamples,
    optimizeRates,
    scenarios);
  for (auto &scenario: scenarios) {
    scenario.gatherReconciliationStatistics(events);
  }
  ParallelContext::barrier();
  events.parallelSum();
}
  

void Routines::getTransfersFrequencies(const std::string &speciesTreeFile,
    Families &families,
    const ModelParameters &modelParameters,
    unsigned int reconciliationSamples,
    TransferFrequencies &transferFrequencies,
    const std::string &outputDir)
{
  const bool bestReconciliation = (reconciliationSamples == 0);
  const bool saveTransfersOnly = true; 
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
  inferReconciliation(speciesTreeFile, 
      families, 
      transfersModelParameter, 
      outputDir, bestReconciliation, 
      reconciliationSamples, 
      optimizeRates, 
      saveTransfersOnly);
  
  SpeciesTree speciesTree(speciesTreeFile);
  for (int i = (bestReconciliation ? -1 : 0); i < (int)reconciliationSamples; ++i) {
    auto begin = ParallelContext::getBegin(families.size());
    auto end = ParallelContext::getEnd(families.size());
    for (unsigned int j = begin; j < end; ++j) {
      auto &family = families[j];
      std::string transfersFile = getTransfersFile(outputDir, family.name, i);
      std::string line;
      std::ifstream is(transfersFile);
      if (!is.good()) {
        Logger::info << "ERROR " << transfersFile << std::endl;
      }
      assert(is.good());
      while (std::getline(is, line))
      {
        std::istringstream iss(line);
        std::string label1;
        std::string label2;
        iss >> label1 >> label2;
        std::string key = getTransferKey(label1, label2);
        if (transferFrequencies.end() == transferFrequencies.find(key)) {
          transferFrequencies.insert(std::pair<std::string, unsigned int>(key, 1));
        } else {
          transferFrequencies[key]++;
        }
      }
    }
  }
  mpiMergeTransferFrequencies(transferFrequencies, outputDir);
  ParallelContext::barrier();
  assert(ParallelContext::isRandConsistent());
}


void Routines::getParametersFromTransferFrequencies(const std::string &speciesTreeFile,
  const TransferFrequencies &frequencies, 
  Parameters &parameters)
{
  SpeciesTree speciesTree(speciesTreeFile);
  std::unordered_map<std::string, unsigned int> labelsToIds;
  speciesTree.getLabelsToId(labelsToIds);
  unsigned int species = speciesTree.getTree().getNodesNumber();
  parameters = Parameters(species * species);
  for (auto &frequency: frequencies) {
    std::string label1;
    std::string label2;
    getLabelsFromTransferKey(frequency.first, label1, label2);
    unsigned int id1 = labelsToIds[label1];
    unsigned int id2 = labelsToIds[label2];
    parameters[id1 * species + id2] = frequency.second;
  }
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

