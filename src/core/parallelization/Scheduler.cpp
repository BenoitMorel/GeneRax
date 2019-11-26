
#include "Scheduler.hpp"
#include <../../ext/MPIScheduler/src/mpischeduler.hpp>
#include <vector>
#include <parallelization/ParallelContext.hpp>
#include <IO/FileSystem.hpp>
#include <cassert>

void Scheduler::schedule(const std::string &outputDir, 
    const std::string &commandFile, 
    bool splitImplem, 
    const std::string &execPath)
{
  assert(ParallelContext::isRandConsistent());
  auto consistentSeed = rand();
  std::vector<char *> argv;
  std::string exec = "mpi-scheduler";
  std::string implem = splitImplem ? "--split-scheduler" : "--fork-scheduler";
  
  std::string called_library = splitImplem ? "--static_scheduled_main" :  execPath;
  std::string jobFailureFatal = "1";
  std::string threadsArg;
  std::string outputLogs = FileSystem::joinPaths(outputDir, "logs.txt");
  std::string ranks = std::to_string(ParallelContext::getSize());
  argv.push_back(const_cast<char *>(exec.c_str()));
  argv.push_back(const_cast<char *>(implem.c_str()));
  argv.push_back(const_cast<char *>(ranks.c_str()));
  argv.push_back(const_cast<char *>(called_library.c_str()));
  argv.push_back(const_cast<char *>(commandFile.c_str()));
  argv.push_back(const_cast<char *>(outputDir.c_str()));
  argv.push_back(const_cast<char *>("--jobs-failure-fatal"));
  argv.push_back(const_cast<char *>("--logs"));
  argv.push_back(const_cast<char *>(outputLogs.c_str()));
  ParallelContext::barrier(); 
  if (splitImplem || ParallelContext::getRank() == 0) {
    void *comm = splitImplem ? static_cast<void *>(&ParallelContext::getComm()) : 0;
    // this is a terrible hack
    std::cout.setstate(std::ios::failbit);
    mpi_scheduler_main(static_cast<int>(argv.size()), &argv[0], comm);
    std::cout.clear();
  }
  // seed accross ranks might not be consistent anymore
  srand(consistentSeed);
  ParallelContext::barrier(); 
}

