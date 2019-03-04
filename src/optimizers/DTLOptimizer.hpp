#pragma once

class JointTree;

class DTLOptimizer {
public:
  static void optimizeDLRates(JointTree &jointTree);
  static void optimizeDTLRates(JointTree &jointTree);
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
  static void optimizeDLRateMN(JointTree &jointTree);
};


