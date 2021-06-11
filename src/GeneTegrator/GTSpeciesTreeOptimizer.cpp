#include "GTSpeciesTreeOptimizer.hpp" 
#include <IO/Logger.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <util/Paths.hpp>

double GTSpeciesTreeLikelihoodEvaluator::computeLikelihood()
{
  double sumLL = 0.0;
  for (auto &evaluation: _evaluations) {
    auto ll = evaluation->computeLogLikelihood();
    sumLL += ll;
  }
  ParallelContext::sumDouble(sumLL);
  return sumLL;
}

void GTSpeciesTreeLikelihoodEvaluator::fillPerFamilyLikelihoods(
    PerFamLL &perFamLL)
{
  ParallelContext::barrier();
  std::vector<double> localLL;
  for (auto &evaluation: _evaluations) {
    localLL.push_back(evaluation->computeLogLikelihood());
  }
  ParallelContext::concatenateHetherogeneousDoubleVectors(
      localLL, perFamLL);
}

GTSpeciesTreeOptimizer::GTSpeciesTreeOptimizer(
    const std::string speciesTreeFile, 
    const Families &families, 
    const RecModelInfo &info,
    const std::string &outputDir):
  _speciesTree(std::make_unique<SpeciesTree>(speciesTreeFile)),
  _outputDir(outputDir),
  _bestRecLL(-std::numeric_limits<double>::infinity())
{
  saveCurrentSpeciesTreeId("starting_species_tree.newick");
  saveCurrentSpeciesTreeId();
  PerCoreGeneTrees perCoreGeneTrees(families, false);
  Logger::timed << "Initializing ccps" << std::endl;
  for (const auto &geneTree: perCoreGeneTrees.getTrees()) {
    auto &family = families[geneTree.familyIndex];
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    switch (info.model) {
    case RecModel::UndatedDL:
      _evaluations.push_back(
        std::make_shared<UndatedDLMultiModel<ScaledValue> >(
        _speciesTree->getTree(),
        mapping,
        info,
        family.startingGeneTree));
      break;
    case RecModel::UndatedDTL:
      _evaluations.push_back(std::make_shared<UndatedDTLMultiModel<ScaledValue> >(
        _speciesTree->getTree(),
        mapping,
        info,
        family.startingGeneTree));
      break;
    default:
      assert(false);
      break;
    }
  }
  Logger::timed << "Initializing ccps finished" << std::endl;
  _speciesTree->addListener(this);
  ParallelContext::barrier();
}


double GTSpeciesTreeOptimizer::computeRecLikelihood()
{
  double sumLL = 0.0;
  for (auto &evaluation: _evaluations) {
    auto ll = evaluation->computeLogLikelihood();
    sumLL += ll;
  }
  ParallelContext::sumDouble(sumLL);
  return sumLL;
}

double GTSpeciesTreeOptimizer::sprSearch(unsigned int radius)
{
  double bestLL = computeRecLikelihood();
  Logger::info << std::endl;
  Logger::timed << "[Species search]" << " Starting species tree local SPR search, radius=" 
    << radius << ", bestLL=" << bestLL  << ", hash=" << _speciesTree->getHash() << ")" <<std::endl;
  double newLL = bestLL;
  do {
    bestLL = newLL;
    newLL = fastSPRRound(radius);
  } while (newLL - bestLL > 0.001);
  bestLL = newLL;
  Logger::timed << "After normal search: LL=" << bestLL << std::endl;
  return bestLL;
}

double GTSpeciesTreeOptimizer::fastSPRRound(unsigned int radius)
{
  Logger::timed << "[Species search] Start SPR Round radius=" << radius << std::endl;
  _bestRecLL = computeRecLikelihood();
  auto hash1 = _speciesTree->getNodeIndexHash(); 
  auto supportValues = std::vector<double>();//_getSupport();
  double maxSupport = 0.2; // ignored for now
  std::vector<unsigned int> prunes;
  SpeciesTreeOperator::getPossiblePrunes(*_speciesTree, 
      prunes,
      supportValues,
      maxSupport);
  for (auto prune: prunes) {
    std::vector<unsigned int> regrafts;
    SpeciesTreeOperator::getPossibleRegrafts(*_speciesTree, prune, radius, regrafts);
    for (auto regraft: regrafts) {
      if (testPruning(prune, regraft)) {
        Logger::timed << "\tbetter tree (LL=" 
          << _bestRecLL << ", hash=" 
          << _speciesTree->getHash() 
          << ") "
          << std::endl;
        hash1 = _speciesTree->getNodeIndexHash(); 
        assert(ParallelContext::isIntEqual(hash1));
        //_bestRecLL = veryLocalSearch(prune);
      }
    }
  }
  return _bestRecLL;
}

void GTSpeciesTreeOptimizer::onSpeciesTreeChange(const std::unordered_set<pll_rnode_t *> *nodesToInvalidate)
{
  for (auto &evaluation: _evaluations) {
    evaluation->onSpeciesTreeChange(nodesToInvalidate);
  }
}

bool GTSpeciesTreeOptimizer::testPruning(unsigned int prune,
    unsigned int regraft)
{
  // Apply the move
  auto rollback = SpeciesTreeOperator::applySPRMove(*_speciesTree, prune, regraft);
  double newLL = computeRecLikelihood();
  if (newLL > _bestRecLL) {
    // Better tree found! keep it and return
    // without rollbacking
    newBestTreeCallback(newLL);
    return true;
  }
  SpeciesTreeOperator::reverseSPRMove(*_speciesTree, prune, rollback);
  return false;
}

void GTSpeciesTreeOptimizer::newBestTreeCallback(double newLL)
{
  saveCurrentSpeciesTreeId();
  _bestRecLL = newLL;
}

std::string GTSpeciesTreeOptimizer::saveCurrentSpeciesTreeId(std::string name, bool masterRankOnly)
{
  std::string res = Paths::getSpeciesTreeFile(_outputDir, name);
  
  saveCurrentSpeciesTreePath(res, masterRankOnly);
  return res;
}

void GTSpeciesTreeOptimizer::saveCurrentSpeciesTreePath(const std::string &str, bool masterRankOnly)
{
  _speciesTree->saveToFile(str, masterRankOnly);
  if (masterRankOnly) {
    ParallelContext::barrier();
  }
}
  
double GTSpeciesTreeOptimizer::rootSearch(unsigned int maxDepth)
{
  GTSpeciesTreeLikelihoodEvaluator evaluator(_evaluations);
  SpeciesRootSearch::rootSearch(
      *_speciesTree,
      evaluator,
      maxDepth);
  saveCurrentSpeciesTreeId();
  return computeRecLikelihood();
}

