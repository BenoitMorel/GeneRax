#include "RaxmlMaster.hpp"

#include <IO/Logger.hpp>
#include <IO/FileSystem.hpp>
#include <IO/LibpllParsers.hpp>
#include <IO/ParallelOfstream.hpp>
#include <parallelization/Scheduler.hpp>
#include <sstream>

void RaxmlMaster::runRaxmlOptimization(Families &families,
    const std::string &output,
    const std::string &execPath,
    int iteration,
    bool splitImplem,
    long &sumElapsed)

{
  auto start = Logger::getElapsedSec();
  Logger::timed << "Starting raxml light step" << std::endl;
  std::stringstream outputDirName;
  outputDirName << "raxml_light_" << iteration;
  std::string outputDir = FileSystem::joinPaths(output, outputDirName.str());
  FileSystem::mkdir(outputDir, true);
  std::string commandFile = FileSystem::joinPaths(outputDir, "raxml_light_command.txt");
  auto geneTreeSizes = LibpllParsers::parallelGetTreeSizes(families);
  ParallelOfstream os(commandFile);
  for (size_t i = 0; i < families.size(); ++i) {
    auto &family = families[i];
    std::string familyOutput = FileSystem::joinPaths(output, "results");
    familyOutput = FileSystem::joinPaths(familyOutput, family.name);
    std::string geneTreePath = FileSystem::joinPaths(familyOutput, "geneTree.newick");
    std::string libpllModelPath = FileSystem::joinPaths(familyOutput, "libpllModel.txt");
    std::string outputStats = FileSystem::joinPaths(familyOutput, "raxml_light_stats.txt");
    auto taxa = geneTreeSizes[i];
    os << family.name << " ";
    os << 1 << " "; // cores
    os << taxa << " " ; // cost
    os << "raxmlLight" << " ";
    os << family.startingGeneTree << " ";
    os << family.alignmentFile << " ";
    os << family.libpllModel  << " ";
    os << geneTreePath << " ";
    os << libpllModelPath << " ";
    os << outputStats <<  std::endl;
    family.startingGeneTree = geneTreePath;
    family.statsFile = outputStats;
    family.libpllModel = libpllModelPath;
  }    
  os.close();
  Scheduler::schedule(outputDir, commandFile, splitImplem, execPath); 
  auto elapsed = (Logger::getElapsedSec() - start);
  sumElapsed += elapsed;
  Logger::timed << "End of raxml light step (after " << elapsed << "s)"  << std::endl;
}

