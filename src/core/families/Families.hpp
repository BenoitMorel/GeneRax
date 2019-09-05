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
  }
};

typedef std::vector<FamilyInfo> Families;
void filterFamilies(Families &families, const std::string &speciesTreeFile);


