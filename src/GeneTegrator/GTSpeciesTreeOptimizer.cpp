#include "GTSpeciesTreeOptimizer.hpp" 
#include <IO/Logger.hpp>


GTSpeciesTreeOptimizer::GTSpeciesTreeOptimizer(
    const std::string speciesTreeFile, 
    const Families &families, 
    const std::string &outputDir):
  _speciesTree(speciesTreeFile)
{
  for (auto &family: families) {
    GeneSpeciesMapping mapping;
    mapping.fill(family.mappingFile, family.startingGeneTree);
    _evaluations.push_back(std::make_shared<MultiEvaluation>(
        _speciesTree,
        mapping,
        family.startingGeneTree));
  }
}
  
double GTSpeciesTreeOptimizer::computeLogLikelihood()
{
  double sumLL = 0.0;
  for (auto &evaluation: _evaluations) {
    auto ll = evaluation->computeLogLikelihood();
    Logger::info << "ll=" << ll << std::endl;
    sumLL += ll;
  }
  return sumLL;
  Logger::info << "total ll = " << sumLL << std::endl;
}

