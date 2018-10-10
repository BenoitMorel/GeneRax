#ifndef _SPR_SEARCH_H_
#define _SPR_SEARCH_H_

class JointTree;

class SPRSearch {
public:
    static void applySPRSearch(JointTree &jointTree);
    static bool applySPRRound(JointTree &jointTree, int radius, double &bestLoglk);
};

#endif
