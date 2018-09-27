#ifndef _NNI_SEARCH_H_
#define _NNI_SEARCH_H_

class ParallelJointTree;

class NNISearch {
public:
    static void applyNNISearch(ParallelJointTree &jointTree);
    static bool applyNNIRound(ParallelJointTree &jointTree, double &bestLoglk);
};

#endif
