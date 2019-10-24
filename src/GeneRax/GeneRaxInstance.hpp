#pragma once


#include "GeneRaxArguments.hpp"
#include <string>
#include <maths/Parameters.hpp>
#include <IO/FamiliesFileParser.hpp>


struct GeneRaxInstance {
  GeneRaxArguments args;
  std::string speciesTree;
  Families initialFamilies;
  Families currentFamilies;
  RecModel recModel;
  Parameters rates;
  double totalLibpllLL;
  double totalRecLL;
  long elapsedRates;
  long elapsedSPR;
  long elapsedRaxml;
  unsigned int currentIteration;
  
  GeneRaxInstance(int argc, char** argv):
    args(argc, argv),
    recModel(ArgumentsHelper::strToRecModel(args.reconciliationModelStr)),
    totalLibpllLL(0),
    totalRecLL(0),
    elapsedRates(0),
    elapsedSPR(0),
    elapsedRaxml(0),
    currentIteration(0)
  {
    switch (recModel) {
      case RecModel::UndatedDL:
      rates = Parameters(2);
      rates[0] = args.dupRate;
      rates[1] = args.lossRate;
      break;
      case RecModel::UndatedDTL: case RecModel::UndatedDTLAdvanced:
      rates = Parameters(3);
      rates[0] = args.dupRate;
      rates[1] = args.lossRate;
      rates[2] = args.transferRate;
      break; 
    }
  }

};


