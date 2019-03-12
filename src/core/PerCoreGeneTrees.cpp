#include "PerCoreGeneTrees.hpp"
#include <ParallelContext.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/Logger.hpp>

PerCoreGeneTrees::PerCoreGeneTrees(const vector<FamiliesFileParser::FamilyInfo> &families)
{
  Logger::error << "We should move some libpll stuff to another file" << endl;
  for (unsigned int i = ParallelContext::getRank(); i < families.size(); i += ParallelContext::getSize()) {
    _geneTrees.push_back(GeneTree());
    _geneTrees.back().mapping = GeneSpeciesMapping(families[i].mappingFile);
    _geneTrees.back().tree = LibpllEvaluation::readNewickFromFile(families[i].startingGeneTree);
  }
}
