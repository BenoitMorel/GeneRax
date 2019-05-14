#pragma once

#include <string>
#include <vector>




class FamiliesFileParser {
public:
  struct FamilyInfo {
    std::string name;
    std::string startingGeneTree;
    std::string alignmentFile;
    std::string mappingFile;
    std::string libpllModel;
    std::string statsFile;
    FamilyInfo() {
      reset();
    }
    void reset() {
      name = "";
      startingGeneTree = "";
      alignmentFile = "";
      mappingFile = "";
      libpllModel = "GTR";
      statsFile = "";
    }
  };
  static std::vector<FamilyInfo> parseFamiliesFile(const std::string &familiesFile);
};
