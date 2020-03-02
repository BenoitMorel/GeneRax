#pragma once

#include <string>


/**
 *  Schedule jobs with the external dependency MPIScheduler, 
 *  allowing to schedule independant jobs in parallel with a 
 *  specified number of cores per job. 
 *  @param outputDir the GeneRax run output directory
 *  @param commandFile command file path (see MPIScheduler 
 *                     documentation)
 *  @param splitImplem split or fork implementaiton (see 
 *                     MPIScheduler documentation)
 *  @param execPath the path to the executable to schedule
 *                  (in practice, it will be the same as the
 *                  main executable, but with specific arguments)
 */
class Scheduler {
public:
  Scheduler() = delete;
  static void schedule(const std::string &outputDir, 
      const std::string &commandFile, 
      bool splitImplem, 
      const std::string &execPath);
};
