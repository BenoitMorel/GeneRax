#pragma once

#include <string>
#include <vector>
using namespace std;



class FamiliesFileParser {
public:
  struct FamilyInfo {
    string name;
    string startingGeneTree;
    string alignmentFile;
    string mappingFile;
    string libpllModel;
    string statsFile;
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
  static vector<FamilyInfo> parseFamiliesFile(const string &familiesFile);
};
