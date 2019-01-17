#pragma once

#include <string>
using namespace std;


class Arguments {
public:
  enum ReconciliationModel {
    UndatedDL,
    UndatedDTL,
    DatedDL,
    InvalidModel
  };
  
public:
  static void init(int argc, char * argv[]);
  static void checkInputs();
  static void printHelp();
  static void printCommand();
  static void printSummary();
public:
  static int argc;
  static char ** argv;
  static string geneTree;
  static string alignment;
  static string speciesTree;
  static string geneSpeciesMap;
  static string strategy;
  static ReconciliationModel reconciliationModel;
  static string output;
  static bool check;
  static bool rootedGeneTree;
  static bool noFelsensteinLikelihood;
  static bool userDTLRates;
  static double dupRate;
  static double lossRate;
  static double transferRate;
};

