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
    color = -1;
  }
};

typedef std::vector<FamilyInfo> Families;
void filterFamilies(Families &families, const std::string &speciesTreeFile);
void duplicatesFamilies(const Families &families, Families &duplicatedFamilies, unsigned int factor);
void contractFamilies(const Families &duplicatedFamilies, Families &families);

void splitInitialFamilies(const Families &families, std::vector<Families> &splitFamilies, unsigned int duplicates, unsigned int splitsNumber);
void mergeSplitFamilies(const std::vector<Families> &splitFamilies, Families &families, unsigned int splitsNumber);


