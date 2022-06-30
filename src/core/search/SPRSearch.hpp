#pragma once

class JointTree;

class SPRSearch {
public:
  virtual ~SPRSearch() {}
    static void applySPRSearch(JointTree &jointTree);
    static bool applySPRRound(JointTree &jointTree, int radius, double &bestLoglk, bool blo = true);
    static bool applySPRRoundDigg(JointTree &jointTree, int radius, double &bestLoglk, bool blo = true);
};

