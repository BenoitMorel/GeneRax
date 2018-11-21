
#ifndef JOINTSEARCH_RECONCILIATIONEVALUATION_HPP
#define JOINTSEARCH_RECONCILIATIONEVALUATION_HPP

#include <Arguments.hpp>
#include <parsers/GeneSpeciesMapping.hpp>
#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <memory>

using namespace std;


class ReconciliationEvaluation {

public:
  ReconciliationEvaluation(pll_rtree_t *speciesTree,
    const GeneSpeciesMapping& map,
    Arguments::ReconciliationModel model);

  void setRates(double dupRate, double lossRate, 
    double transfers = 0.0);

  pll_unode_t * getRoot() {return reconciliationModel->getRoot();}
  void setRoot(pll_unode_t *root) {reconciliationModel->setRoot(root);}
  double evaluate(shared_ptr<pllmod_treeinfo_t> treeinfo);

private:
  shared_ptr<AbstractReconciliationModel> reconciliationModel;
  bool firstCall;
};

#endif
