#pragma once

class JointTree;
class PerCoreGeneTrees;
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <string>
#include <enums.hpp>

using namespace std;


class DTLOptimizer {
public:
  static void optimizeDLRates(JointTree &jointTree, RecOpt method);
  static void optimizeDTLRates(JointTree &jointTree, RecOpt method);
  
  static void optimizeDLRates(PerCoreGeneTrees &geneTrees, pll_rtree_t *speciesTree, RecModel model);
private:
  static void findBestRatesDTL(JointTree &jointTree,
      double minDup, double maxDup,
      double minLoss, double maxLoss, 
      double minTrans, double maxTrans, 
      int steps,
      double &bestDup,
      double &bestLoss,
      double &bestTrans,
      double &bestLL);

  static void findBestRatesDL(JointTree &jointTree, 
      double minDup, double maxDup,
      double minLoss, double maxLoss, int steps,
      double &bestDup,
      double &bestLoss,
      double &bestLL) ;
  static void optimizeRateSimplex(JointTree &jointTree, bool transfers);
  static void optimizeDLRatesWindow(JointTree &jointTree);
  static void optimizeDTLRatesWindow(JointTree &jointTree);
};


