#ifndef _NNI_SEARCH_H_
#define _NNI_SEARCH_H_


#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <ale/containers/GeneMap.h>
#include <treeSearch/JointTree.h>
#include <treeSearch/Moves.h>
#include <ale/tools/SpeciesGeneMapper.h>



class NNISearch {
public:
    static void applyNNISearch(AbstractJointTree &jointTree);
    static bool applyNNIRound(AbstractJointTree &jointTree, double &bestLoglk);
};

#endif
