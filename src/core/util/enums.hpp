#pragma once

#include <cassert>


/**
 *  Reconciliation models 
 */
enum class RecModel {
  UndatedDL, UndatedDTL
};


/*
 *  DTLRates numerical optimization methods
 */
enum class RecOpt {
  Grid, Simplex, Gradient
};


/*
 * GeneRax search modes
 */
enum class Strategy {
  SPR, EVAL
};


enum class SpeciesStrategy {
  SPR, TRANSFERS, HYBRID
};

/*
 *  Output formats for reconciled gene trees
 */
enum class ReconciliationFormat {
  NHX = 0, RecPhyloXML
};


enum class ReconciliationEventType {
  EVENT_S = 0 , EVENT_SL, EVENT_D, EVENT_T, EVENT_TL, EVENT_L, EVENT_None, EVENT_Invalid
};


/*
 * Defines how to reuse computations when computing
 * the reconciliation likelihood
 */
enum class PartialLikelihoodMode {
  PartialGenes = 0, // reuse per-gene CLVs 
  PartialSpecies, // reuse per-species CLVs
  NoPartial // always recompute all CLVs from scratch
};

class Enums {
public:
  Enums() = delete;

  static unsigned int freeParameters(RecModel m)  {
    switch (m) {
      case RecModel::UndatedDL:
        return 2;
      case RecModel::UndatedDTL:
        return 3;
    }
    assert(false);
  }

  static bool accountsForTransfers(RecModel m) 
  {
    switch (m) {
    case RecModel::UndatedDL:
      return false;
    case RecModel::UndatedDTL:
      return true;
    }
    assert(false);
    return false;
  }

  static bool implementsApproxLikelihood(RecModel m)
  {
    switch (m) {
      case RecModel::UndatedDL:
        return false;
      case RecModel::UndatedDTL:
        return true;
    }
    assert(false);
    return false;
  }

};


