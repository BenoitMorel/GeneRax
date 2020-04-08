#include "GeneRaxInstance.hpp"
#include <IO/FileSystem.hpp>

void GeneRaxInstance::readModelParameters(ModelParameters &modelParameters)
{
  Parameters defaultRates = rates;
  modelParameters = ModelParameters(defaultRates, 
      recModel, 
      args.perFamilyDTLRates, 
      currentFamilies.size());
  auto freeParameters = Enums::freeParameters(recModel);
  Parameters ratesToRead(rates);
  auto resultsPath = FileSystem::joinPaths(args.output, "results");
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
}


