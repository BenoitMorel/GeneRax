#pragma once

class JointTree;

class SPRSearch {
public:
    static void applySPRSearch(JointTree &jointTree);
    static bool applySPRRound(JointTree &jointTree, int radius, double &bestLoglk);
};

