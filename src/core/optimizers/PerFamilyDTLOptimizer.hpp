#pragma once

class JointTree;
#include <likelihoods/ReconciliationEvaluation.hpp>
#include <string>
#include <util/enums.hpp>
#include <maths/DTLRates.hpp>




class PerFamilyDTLOptimizer {
public:
  /**
   * Per-family rates optimization
   */
  static void optimizeDLRates(JointTree &jointTree, RecOpt method);
  static void optimizeDTLRates(JointTree &jointTree, RecOpt method);
 
private:
  static void findBestRatesDTL(JointTree &jointTree,
      double minDup, double maxDup,
      double minLoss, double maxLoss, 
      double minTrans, double maxTrans, 
      unsigned int steps,
      double &bestDup,
      double &bestLoss,
      double &bestTrans,
      double &bestLL);

  static void findBestRatesDL(JointTree &jointTree, 
      double minDup, double maxDup,
      double minLoss, double maxLoss, unsigned int steps,
      double &bestDup,
      double &bestLoss,
      double &bestLL) ;
  static void optimizePerSpeciesRateSimplex(JointTree &jointTree, bool transfers);
  static void optimizeRateSimplex(JointTree &jointTree, bool transfers);
  static void optimizeDLRatesWindow(JointTree &jointTree);
  static void optimizeDTLRatesWindow(JointTree &jointTree);
};



