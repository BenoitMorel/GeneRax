#include "GeneRaxMaster.hpp"

#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/LibpllParsers.hpp>
#include <IO/ParallelOfstream.hpp>
#include <maths/Parameters.hpp>
#include <parallelization/Scheduler.hpp>
#include <sstream>


std::string toArg(const std::string &str) {
  return str.size() ? str : "NONE";
}


void GeneRaxMaster::optimizeGeneTrees(Families &families,
    RecModel recModel,
    Parameters &rates,
    const std::string &output,
    const std::string &resultName,
    const std::string &execPath, 
    const std::string &speciesTreePath,
    RecOpt recOpt,
    bool perFamilyDTLRates,
    bool rootedGeneTree,
    double supportThreshold,
    double recWeight,
    bool enableRec,
    bool enableLibpll,
    int sprRadius,
    int iteration,
    bool schedulerSplitImplem,
    long &elapsed,
    bool inPlace) 
{
  auto start = Logger::getElapsedSec();
  std::stringstream outputDirName;
  outputDirName << "gene_optimization_" << iteration;
  std::string outputDir = FileSystem::joinPaths(output, outputDirName.str());
  FileSystem::mkdir(outputDir, true);
  std::string commandFile = FileSystem::joinPaths(outputDir, "opt_genes_command.txt");
  auto geneTreeSizes = LibpllParsers::parallelGetTreeSizes(families);
  ParallelOfstream os(commandFile);
  std::string ratesFile = FileSystem::joinPaths(outputDir, "dtl_rates.txt");
  rates.save(ratesFile);
  for (size_t i = 0; i < families.size(); ++i) {
    auto &family = families[i];
    std::string familyOutput = FileSystem::joinPaths(output, resultName);
    familyOutput = FileSystem::joinPaths(familyOutput, family.name);
    std::string geneTreePath = FileSystem::joinPaths(familyOutput, "geneTree.newick");
    if (inPlace) {
      // todobenoit make this the normal behavior?
      geneTreePath = family.startingGeneTree;
    }
    std::string outputStats = FileSystem::joinPaths(familyOutput, "stats.txt");
    auto taxa = geneTreeSizes[i];
    unsigned int cores = 1;
    if (sprRadius == 1) {
      cores = taxa / 2;;
    } else if (sprRadius == 2) {
      cores = taxa / 2;
    } else if (sprRadius >= 3) {
      cores = taxa;
    }
    if (cores == 0) {
      cores = 1;
    }
    if (!schedulerSplitImplem) {
      cores = 1;
    }
    os << family.name << " ";
    os << cores << " "; // cores
    os << taxa << " " ; // cost
    os << "optimizeGeneTrees" << " ";
    os << family.startingGeneTree << " ";
    os << toArg(family.mappingFile) << " ";
    os << family.alignmentFile << " ";
    os << speciesTreePath << " ";
    os << family.libpllModel  << " ";
    os << ratesFile << " ";
    os << static_cast<int>(recModel)  << " ";
    os << static_cast<int>(recOpt)  << " ";
    os << static_cast<int>(perFamilyDTLRates)  << " ";
    os << static_cast<int>(rootedGeneTree)  << " ";
    os << supportThreshold << " ";
    os << recWeight  << " ";
    os << static_cast<int>(enableRec)  << " ";
    os << static_cast<int>(enableLibpll)  << " ";
    os << sprRadius  << " ";
    os << geneTreePath << " ";
    os << outputStats <<  std::endl;
    family.startingGeneTree = geneTreePath;
    family.statsFile = outputStats;
  } 
  os.close();
  Scheduler::schedule(outputDir, commandFile, schedulerSplitImplem, execPath); 
  elapsed = (Logger::getElapsedSec() - start);
}

