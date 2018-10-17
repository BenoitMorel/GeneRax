
#ifndef JOINTSEARCH_ALEEVALUATION_HPP
#define JOINTSEARCH_ALEEVALUATION_HPP

// Include ALE
#include <ale/tools/ALE/ALE.h>
#ifdef DTL
#include <ale/tools/ALE/exODT.h>
#else
#include <ale/tools/ALE/exODT_DL.h>
#endif
#include <ale/tools/ALE/fractionMissing.h>
#include <ale/tools/SpeciesGeneMapper.h>
#include <ale/tools/IO/IO.h>

class ALEEvaluation {

public:
  ALEEvaluation(const bpp::PhyloTree& speciestree,
    const SpeciesGeneMap& map);

  void setRates(double dupRate, double lossRate);

  double evaluate(const bpp::PhyloTree& genetree);
private:
#ifdef DTL
  exODT_model model;
#else
  exODT_DL_model model;
#endif
};


#endif //TREERECS_ALEEVALUATION_H
