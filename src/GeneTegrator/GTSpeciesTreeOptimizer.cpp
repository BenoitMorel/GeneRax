#include "GTSpeciesTreeOptimizer.hpp" 
#include <IO/Logger.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>

GTSpeciesTreeOptimizer::GTSpeciesTreeOptimizer(
    const std::string speciesTreeFile, 
    const Families &families, 
    const std::string &outputDir):
  _speciesTree(speciesTreeFile)
{
  PerCoreGeneTrees perCoreGeneTrees(families, false);
  for (const auto &geneTree: perCoreGeneTrees.getTrees()) {
    auto &family = families[geneTree.familyIndex];
  //for (auto &family: families) {
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    _evaluations.push_back(std::make_shared<MultiEvaluation>(
        _speciesTree,
        mapping,
        family.startingGeneTree));
  }

  ParallelContext::barrier();
}
  
double GTSpeciesTreeOptimizer::computeLogLikelihood()
{
  double sumLL = 0.0;
  for (auto &evaluation: _evaluations) {
    auto ll = evaluation->computeLogLikelihood();
    Logger::info << "ll=" << ll << std::endl;
    sumLL += ll;
  }
  ParallelContext::sumDouble(sumLL);
  return sumLL;
}

