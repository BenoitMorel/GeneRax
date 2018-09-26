#ifndef _SPR_SEARCH_H_
#define _SPR_SEARCH_H_


#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <ale/containers/GeneMap.h>
#include <treeSearch/JointTree.h>
#include <treeSearch/Moves.h>

class SPRSearch {
public:
    static void applySPRSearch(AbstractJointTree &jointTree);
    static bool applySPRRound(AbstractJointTree &jointTree, int radius, double &bestLoglk);

};

#endif
