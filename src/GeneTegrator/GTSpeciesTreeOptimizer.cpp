#include "GTSpeciesTreeOptimizer.hpp" 
#include <IO/Logger.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <util/Paths.hpp>
#include <search/SpeciesSPRSearch.hpp>

double GTSpeciesTreeLikelihoodEvaluator::computeLikelihood()
{
  return computeLikelihoodFast();
}

double GTSpeciesTreeLikelihoodEvaluator::computeLikelihoodFast()
{
  double sumLL = 0.0;
  for (auto &evaluation: *_evaluations) {
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
  for (auto &evaluation: *_evaluations) {
    localLL.push_back(evaluation->computeLogLikelihood());
  }
  ParallelContext::concatenateHetherogeneousDoubleVectors(
      localLL, perFamLL);
}

void GTSpeciesTreeLikelihoodEvaluator::countEvents()
{
  std::vector<unsigned int> events(
      static_cast<unsigned int>(ReconciliationEventType::EVENT_Invalid),
      0);
  for (auto &evaluation: *_evaluations) {
    Scenario scenario;
    bool stochastic = true;
    evaluation->inferMLScenario(scenario, stochastic);
    for (unsigned int i = 0; i < events.size(); ++i) {
      events[i] += scenario.getEventCount(
          static_cast<ReconciliationEventType>(i));
    }
  }
  ParallelContext::sumVectorUInt(events);
  Logger::info << "S=" << events[(unsigned int)ReconciliationEventType::EVENT_S] << std::endl;
  Logger::info << "SL=" << events[(unsigned int)ReconciliationEventType::EVENT_SL] << std::endl;
  Logger::info << "D=" << events[(unsigned int)ReconciliationEventType::EVENT_D] << std::endl;
  Logger::info << "T=" << events[(unsigned int)ReconciliationEventType::EVENT_T] << std::endl;
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
  _evaluator.setEvaluations(_evaluations);
  Logger::timed << "Initial ll=" << computeRecLikelihood() << std::endl;
  _evaluator.countEvents();
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
  AverageStream useless;
  if (SpeciesSPRSearch::SPRSearch(*_speciesTree,
      _evaluator,
      useless,
      radius,
      bestLL,
      bestLL)) {
    newBestTreeCallback(bestLL);
  }
  Logger::timed << "After normal search: LL=" << bestLL << std::endl;
  return bestLL;
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
  SpeciesRootSearch::rootSearch(
      *_speciesTree,
      _evaluator,
      maxDepth);
  saveCurrentSpeciesTreeId();
  return computeRecLikelihood();
}

