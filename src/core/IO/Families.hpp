#pragma once

#include <string>
#include <vector>



struct FamilyInfo {
  std::string name;
  std::string startingGeneTree;
  std::string alignmentFile;
  std::string mappingFile;
  std::string libpllModel;
  std::string statsFile;
  unsigned int color;
  FamilyInfo() {
    reset();
  }
  void reset() {
    name = "";
    startingGeneTree = "__random__";
    alignmentFile = "";
    mappingFile = "";
    libpllModel = "GTR";
    statsFile = "";
    color = 0;
  }
};

typedef std::vector<FamilyInfo> Families;

class Family {
public:
  Family() = delete;
  static void filterFamilies(Families &families, const std::string &speciesTreeFile, bool checkAlignments, bool checkSpeciesTree);
  static void printStats(Families &families, const std::string &speciesTreeFile, const std::string &coverageFile, const std::string &fractionMissingFile);
};


