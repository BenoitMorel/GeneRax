
#ifndef JOINTSEARCH_ALEEVALUATION_HPP
#define JOINTSEARCH_ALEEVALUATION_HPP

// Include ALE

#include <ale/tools/SpeciesGeneMapper.h>
#include <parsers/GeneSpeciesMapping.hpp>
#include <likelihoods/ale/UndatedDLModel.hpp>

class ALEEvaluation {

public:
  ALEEvaluation(pll_rtree_t *speciesTree,
    const GeneSpeciesMapping& map);

  void setRates(double dupRate, double lossRate);

  pll_unode_t * getRoot() {return undatedDLModel.getRoot();}
  void setRoot(pll_unode_t *root) {undatedDLModel.setRoot(root);}
  double evaluate(shared_ptr<pllmod_treeinfo_t> treeinfo);
private:
  UndatedDLModel undatedDLModel;
  bool firstCall;
};


#endif //TREERECS_ALEEVALUATION_H