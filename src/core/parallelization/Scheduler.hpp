#pragma once

#include <string>


class Scheduler {
public:
  Scheduler() = delete;
  static void schedule(const std::string &outputDir, const std::string &commandFile, bool splitImplem, const std::string &execPath);
};
