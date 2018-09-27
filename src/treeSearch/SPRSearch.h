#ifndef _SPR_SEARCH_H_
#define _SPR_SEARCH_H_

class ParallelJointTree;

class SPRSearch {
public:
    static void applySPRSearch(ParallelJointTree &jointTree);
    static bool applySPRRound(ParallelJointTree &jointTree, int radius, double &bestLoglk);
};

#endif
