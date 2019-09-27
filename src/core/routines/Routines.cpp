
#include "Routines.hpp"
#include <sstream>
#include <IO/Logger.hpp>
#include <IO/LibpllParsers.hpp>
#include <trees/PerCoreGeneTrees.hpp>
#include <trees/SpeciesTree.hpp>
#include <optimizers/DTLOptimizer.hpp>
#include <parallelization/ParallelContext.hpp>
#include <IO/FileSystem.hpp>
#include <likelihoods/LibpllEvaluation.hpp>

void Routines::optimizeRates(bool userDTLRates, 
    const std::string &speciesTreeFile,
    RecModel recModel,
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
  pll_rtree_t *speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile); 
  if (perSpeciesRates) {
    //auto ratesGlobal = DTLOptimizer::optimizeParametersGlobalDTL(geneTrees, speciesTree, recModel);
    //auto ratesEmpirical = ratesGlobal;
    //optimizeSpeciesRatesEmpirical(speciesTreeFile, recModel, families, ratesEmpirical, outputDir, sumElapsed);
    rates = DTLOptimizer::optimizeParametersPerSpecies(geneTrees, speciesTree, recModel);
    /*
    Logger::timed << " Per-species DTL rates optimization from global rates: RecLL=" << rates1.getScore() << std::endl;
    auto rates2 = DTLOptimizer::optimizeParameters(geneTrees, speciesTree, recModel, ratesEmpirical);
    Logger::timed << " Per-species DTL rates optimization from empitical rates: RecLL=" << rates2.getScore() << std::endl;
    rates = (rates1.getScore() > rates2.getScore() ? rates1 : rates2);
    */
  } else {
    rates = DTLOptimizer::optimizeParametersGlobalDTL(geneTrees, speciesTree, recModel);
  }
  pll_rtree_destroy(speciesTree, 0);
  ParallelContext::barrier(); 
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
}


static std::string getSpeciesEventCountFile(const std::string &outputDir, const std::string &familyName)
{
  return FileSystem::joinPaths(outputDir, FileSystem::joinPaths("reconciliations", familyName + "_speciesEventCounts.txt"));
}

static std::string getTransfersFile(const std::string &outputDir, const std::string &familyName)
{
  return FileSystem::joinPaths(outputDir, FileSystem::joinPaths("reconciliations", familyName + "_transfers.txt"));
}

void Routines::inferReconciliation(
    const std::string &speciesTreeFile,
    Families &families,
    RecModel model,
    const Parameters &rates,
    const std::string &outputDir
    )
{
  ParallelContext::barrier();
  pll_rtree_t *speciesTree = LibpllParsers::readRootedFromFile(speciesTreeFile); 
  PerCoreGeneTrees geneTrees(families);
  std::string reconciliationsDir = FileSystem::joinPaths(outputDir, "reconciliations");
  FileSystem::mkdir(reconciliationsDir, true);
  auto speciesNodesCount = speciesTree->tip_count + speciesTree->inner_count;
  std::vector<double> dup_count(speciesNodesCount, 0.0);
  ParallelContext::barrier();
  for (auto &tree: geneTrees.getTrees()) {
    std::string eventCountsFile = FileSystem::joinPaths(reconciliationsDir, tree.name + "_eventCounts.txt");
    std::string speciesEventCountsFile = getSpeciesEventCountFile(outputDir, tree.name);
    std::string transfersFile = getTransfersFile(outputDir, tree.name);
    std::string treeWithEventsFileNHX = FileSystem::joinPaths(reconciliationsDir, tree.name + "_reconciliated.nhx");
    std::string treeWithEventsFileRecPhyloXML = FileSystem::joinPaths(reconciliationsDir, tree.name + "_reconciliated.xml");
    Scenario scenario;
    ReconciliationEvaluation evaluation(speciesTree, tree.mapping, model, true);
    evaluation.setRates(rates);
    evaluation.evaluate(tree.tree);
    evaluation.inferMLScenario(tree.tree, scenario);
    scenario.saveEventsCounts(eventCountsFile, false);
    scenario.savePerSpeciesEventsCounts(speciesEventCountsFile, false);
    scenario.saveTransfers(transfersFile, false);
    scenario.saveReconciliation(treeWithEventsFileRecPhyloXML, RecPhyloXML, false);
    scenario.saveReconciliation(treeWithEventsFileNHX, NHX, false);
  }
  pll_rtree_destroy(speciesTree, 0);
  ParallelContext::barrier();
}

bool Routines::createRandomTrees(const std::string &geneRaxOutputDir, Families &families)
{
  std::string startingTreesDir = FileSystem::joinPaths(geneRaxOutputDir, "startingTrees");
  bool startingTreesDirCreated = false;
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
  
void Routines::optimizeSpeciesRatesEmpirical(const std::string &speciesTreeFile,
    RecModel recModel,
    Families &families,
    Parameters &rates,
    const std::string &outputDir,
    long &sumElapsed)
{
  ParallelContext::barrier();
  auto start = Logger::getElapsedSec();
  inferReconciliation(speciesTreeFile, families, recModel, rates, outputDir);
  SpeciesTree speciesTree(speciesTreeFile);
  std::unordered_set<std::string> labels;
  std::unordered_map<std::string, std::vector<double> > frequencies;
  speciesTree.getLabels(labels, false);
  std::vector<double> defaultFrequences(4, 0.0);
  for (auto &label: labels) {
    frequencies.insert(std::pair<std::string, std::vector<double>>(label, defaultFrequences));
  }
  for (auto &family: families) {
    std::string speciesEventCountsFile = getSpeciesEventCountFile(outputDir, family.name);
    std::ifstream is(speciesEventCountsFile);
    assert(is.good());
    std::string line;
    while (std::getline(is, line))
    {
      std::istringstream iss(line);
      std::string label;
      iss >> label;
      auto &speciesFreq = frequencies[label];
      for (unsigned int i = 0; i < 4; ++i) {
        double temp;
        iss>> temp;
        speciesFreq[i] += temp;
      }
    }
  }
  unsigned int perSpeciesFreeParameters = Enums::freeParameters(recModel); 
  rates = Parameters(speciesTree.getNodesNumber() * perSpeciesFreeParameters);
  for (unsigned int i = 0; i < speciesTree.getNodesNumber(); ++i) {
    auto speciesNode = speciesTree.getNode(i);
    auto &speciesFreq = frequencies[speciesNode->label];
    double S = speciesFreq[0] + 1.0;
    for (unsigned int j = 0; j < perSpeciesFreeParameters; ++j) {
      unsigned int ratesIndex = speciesNode->node_index * perSpeciesFreeParameters + j;
      rates[ratesIndex] = (speciesFreq[j + 1] + 1.0) / S;
    }
  }
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
  ParallelContext::barrier();
}

static const std::string keyDelimiter("-_-");

std::string getTransferKey(const std::string &label1, const std::string &label2)
{
  return label1 + keyDelimiter + label2; 
}

void getLabelsFromTransferKey(const std::string &key, std::string &label1, std::string &label2)
{
  int pos = key.find(keyDelimiter);
  label1 = key.substr(0, pos);
  label2 = key.substr(pos + keyDelimiter.size());
}
void Routines::getTransfersFrequencies(const std::string &speciesTreeFile,
    RecModel recModel,
    Families &families,
    const Parameters &rates,
    TransferFrequencies &transferFrequencies,
    const std::string &outputDir)
{
  inferReconciliation(speciesTreeFile, families, recModel, rates, outputDir);
  SpeciesTree speciesTree(speciesTreeFile);
  
  for (auto &family: families) {
    std::string transfersFile = getTransfersFile(outputDir, family.name);
    std::string line;
    std::ifstream is(transfersFile);
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
  ParallelOfstream os(FileSystem::joinPaths(outputDir, "transfers.txt"));
  for (auto &freq: transferFrequencies) {
    os << freq.first << " " << freq.second << std::endl;
  }
}


void Routines::getParametersFromTransferFrequencies(const std::string &speciesTreeFile,
  const TransferFrequencies &frequencies, 
  Parameters &parameters)
{
  SpeciesTree speciesTree(speciesTreeFile);
  std::unordered_map<std::string, unsigned int> labelsToIds;
  speciesTree.getLabelsToId(labelsToIds);
  unsigned int species = speciesTree.getNodesNumber();
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








