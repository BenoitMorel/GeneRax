#include "GeneRaxInstance.hpp"
#include <IO/FileSystem.hpp>

void GeneRaxInstance::readModelParameters(ModelParameters &modelParameters)
{
  auto freeParameters = Enums::freeParameters(recModelInfo.model);
  if (!freeParameters) {
    return;
  }
  Parameters defaultRates = rates;
  assert(rates.dimensions());
  
  modelParameters = ModelParameters(defaultRates, 
      currentFamilies.size(),
      getRecModelInfo());
  if (!args.perFamilyDTLRates) {
    return;
  }
  Parameters ratesToRead(rates);
  auto resultsPath = FileSystem::joinPaths(args.outputPath, "results");
  for (unsigned int i = 0; i < this->currentFamilies.size(); ++i) {
    const auto &family = this->currentFamilies[i];
    std::string stats = FileSystem::joinPaths(FileSystem::joinPaths(resultsPath, family.name), "stats.txt");
    std::ifstream paramsReader(stats);
    if (!paramsReader) {
      continue;
    }
    double val = 0.0;
    std::string skipStr;
    // skip  likelihoods
    paramsReader >> val >> val;
    paramsReader >> skipStr >> skipStr >> skipStr;
    for (unsigned int dim = 0; dim < freeParameters; ++dim) {
      paramsReader >> ratesToRead[dim];
    }
    modelParameters.setRates(i, ratesToRead);
  }
  assert(modelParameters.rates.dimensions());
}

RecModelInfo GeneRaxInstance::getRecModelInfo()
{
  return recModelInfo;
}

