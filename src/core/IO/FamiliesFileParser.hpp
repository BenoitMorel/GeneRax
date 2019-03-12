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
    void reset() {
      name = "";
      startingGeneTree = "";
      alignmentFile = "";
      mappingFile = "";
    }
  };
  static vector<FamilyInfo> parseFamiliesFile(const string &familiesFile);
};
