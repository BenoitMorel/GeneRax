#pragma once

#include <util/enums.hpp>

struct RecModelInfo {
  // reconciliation model (UndatedDTL, UndatedDL, etc)
  RecModel model;
  // if set to true, each family can have different set of rates
  bool perFamilyRates;
  // if set to true,  for each family, we prune from the species
  // tree the taxa that are not covered in this family
  bool pruneSpeciesTree;
  bool rootedGeneTree;
  // if the reconciliaiton model accounts for polytomies, branches
  // with a lenghts <= branchLengthThreshold are be contracted
  double branchLengthThreshold;
  
  std::string fractionMissingFile;
  
  
  RecModelInfo():
    model(RecModel::UndatedDTL),
    perFamilyRates(false),
    pruneSpeciesTree(false),
    branchLengthThreshold(-1.0)
  {

  }

  RecModelInfo(RecModel model,
      bool perFamilyRates,
      bool pruneSpeciesTree,
      double branchLengthThreshold,
      const std::string &fractionMissingFile):
    model(model),
    perFamilyRates(perFamilyRates),
    pruneSpeciesTree(pruneSpeciesTree),
    branchLengthThreshold(branchLengthThreshold),
    fractionMissingFile(fractionMissingFile)
  {

  }

  void readFromArgv(char** argv, int &i)
  {
    model = RecModel(atoi(argv[i++]));  
    perFamilyRates = bool(atoi(argv[i++]));
    pruneSpeciesTree = bool(atoi(argv[i++]));
    branchLengthThreshold = double(atof(argv[i++]));
    rootedGeneTree = bool(atoi(argv[i++]));
    fractionMissingFile = std::string(argv[i++]);
    if (fractionMissingFile == "NONE") {
      fractionMissingFile = std::string();
    }
  }

  std::vector<std::string> getArgv() const
  {
    std::vector<std::string> argv;
    argv.push_back(std::to_string(static_cast<int>(model)));
    argv.push_back(std::to_string(static_cast<int>(perFamilyRates)));
    argv.push_back(std::to_string(static_cast<int>(pruneSpeciesTree)));
    argv.push_back(std::to_string(branchLengthThreshold));
    if (fractionMissingFile.size()) {
      argv.push_back(fractionMissingFile);
    } else {
      argv.push_back(std::string("NONE"));
    }
    return argv;
  }

  static int getArgc() 
  {
    return 5;
  }

  unsigned int modelFreeParameters() const {
    return Enums::freeParameters(model);
  }

};
