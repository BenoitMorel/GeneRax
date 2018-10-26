#ifndef _NNI_SEARCH_H_
#define _NNI_SEARCH_H_

class JointTree;

class NNISearch {
public:
    static void applyNNISearch(JointTree &jointTree);
    static bool applyNNIRound(JointTree &jointTree, double &bestLoglk, int bloRadius);
};

#endif
