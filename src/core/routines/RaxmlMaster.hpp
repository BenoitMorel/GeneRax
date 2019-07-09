#pragma once

#include <string>
#include <vector>
#include <IO/FamiliesFileParser.hpp>

class RaxmlMaster {
public:
  static void runRaxmlOptimization(std::vector<FamiliesFileParser::FamilyInfo> &families,
    const std::string &output,
    const std::string &execPath,
    int iteration,
    bool splitImplem,
    long &sumElapsed);

};
