#pragma once

#include <search/SpeciesRootSearch.hpp>
#include <trees/SpeciesTree.hpp>
#include <IO/FamiliesFileParser.hpp>
#include "UndatedDLMultiModel.hpp"
#include "UndatedDTLMultiModel.hpp"
#include <trees/PLLRootedTree.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>
#include <memory>
#include <vector>
#include <maths/ModelParameters.hpp>

class RecModelInfo;
using MultiEvaluation = MultiModel;
using MultiEvaluationPtr = 
  std::shared_ptr<MultiEvaluation>;
using PerCoreMultiEvaluation = std::vector<MultiEvaluationPtr>;

class GTSpeciesTreeLikelihoodEvaluator: public SpeciesTreeLikelihoodEvaluatorInterface {
public:
  GTSpeciesTreeLikelihoodEvaluator(SpeciesTree &speciesTree,
      ModelParameters &modelRates, 
      const Families &families,
      PerCoreGeneTrees &geneTrees);
  virtual ~GTSpeciesTreeLikelihoodEvaluator() {}
  virtual double computeLikelihood();
  virtual double computeLikelihoodFast();
  virtual bool providesFastLikelihoodImpl() const {return false;}
  virtual bool isDated() const {return _modelRates.info.isDated();}
  virtual double optimizeModelRates(bool thorough = false);
  virtual void pushRollback() {}
  virtual void popAndApplyRollback() {}
  virtual void fillPerFamilyLikelihoods(PerFamLL &perFamLL);
  virtual void getTransferInformation(SpeciesTree &speciesTree,
    TransferFrequencies &frequencies,
    PerSpeciesEvents &perSpeciesEvents,
    PerCorePotentialTransfers &potentialTransfers);
  virtual bool pruneSpeciesTree() const {return _modelRates.info.pruneSpeciesTree;}
  virtual void setAlpha(double alpha);
  virtual void onSpeciesDatesChange();  
  virtual void onSpeciesTreeChange(
      const std::unordered_set<corax_rnode_t *> *nodesToInvalidate);
  void printHightPrecisionCount();
  MultiModel &getEvaluation(unsigned int i) {return *_evaluations[i];}


  void addHighway(const Highway &highway);
  void removeHighway();
  void sampleScenarios(unsigned int family, unsigned int samples,
      std::vector<Scenario> &scenarios);
protected:
  virtual double optimizeGammaRates();
  void resetEvaluation(unsigned int i, bool highPrecision);
private:
  SpeciesTree &_speciesTree;
  ModelParameters &_modelRates;
  std::vector<Highway> _highways;
  const Families &_families;
  PerCoreMultiEvaluation _evaluations;
  PerCoreMultiEvaluation _approxEvaluations;
  PerCoreGeneTrees &_geneTrees;
  std::vector<int> _highPrecisions;
  
};


