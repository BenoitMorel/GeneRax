#pragma once


#include "GeneRaxArguments.hpp"
#include <string>
#include <maths/Parameters.hpp>
#include <maths/ModelParameters.hpp>
#include <IO/FamiliesFileParser.hpp>


struct GeneRaxInstance {
  GeneRaxArguments args;
  std::string speciesTree;
  Families initialFamilies;
  Families currentFamilies;
  RecModelInfo recModelInfo;
  Parameters rates;
  ModelParameters modelParameters;
  double totalLibpllLL;
  double totalRecLL;
  long elapsedRates;
  long elapsedSPR;
  long elapsedRaxml;
  unsigned int currentIteration;
  
  GeneRaxInstance(int argc, char** argv):
    args(argc, argv),
    totalLibpllLL(0),
    totalRecLL(0),
    elapsedRates(0),
    elapsedSPR(0),
    elapsedRaxml(0),
    currentIteration(0)
  {
    auto recModel = ArgumentsHelper::strToRecModel(
        args.reconciliationModelStr);
    std::string fractionMissingFile;
    recModelInfo = RecModelInfo(recModel,
        args.perFamilyDTLRates,
        args.gammaCategories,
        args.originationStrategy,
        args.pruneSpeciesTree,
        args.rootedGeneTree,
        args.minGeneBranchLength,
        args.transferConstraint,
        args.noDup,
        args.fractionMissingFile
        );
    rates = getUserParameters();
  }
  
  Parameters getUserParameters() const {
    Parameters user(args.dupRate, args.lossRate, args.transferRate);
    return recModelInfo.getParametersFromUser(user);
  }

  GeneRaxInstance(const GeneRaxInstance &) = delete;
  GeneRaxInstance & operator = (const GeneRaxInstance &) = delete;
  GeneRaxInstance(GeneRaxInstance &&) = delete;
  GeneRaxInstance & operator = (GeneRaxInstance &&) = delete;
  
  void readModelParameters(ModelParameters &modelParameters);
  RecModelInfo getRecModelInfo();
};


