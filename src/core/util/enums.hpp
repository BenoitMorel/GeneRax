#pragma once

#include <cassert>
#include <string>
#include <unordered_map>

using StringToUintMap = std::unordered_map<std::string, unsigned int>;

/**
 *  Reconciliation models 
 */
enum class RecModel {
  UndatedDL, UndatedDTL, UndatedIDTL
};


/*
 *  DTLRates numerical optimization methods
 */
enum class RecOpt {
  Grid, Simplex, Gradient
};


/*
 * Gene tree search mode
 */
enum class GeneSearchStrategy {
  SPR, EVAL
};

/**
 * Species tree search mode
 */
enum class SpeciesSearchStrategy {
  SPR, TRANSFERS, HYBRID
};

/*
 *  Output formats for reconciled gene trees
 */
enum class ReconciliationFormat {
  NHX = 0, RecPhyloXML
};


/**
 * Nature of a reconciliation event
 */ 
enum class ReconciliationEventType {
  EVENT_S = 0,  // speciation
  EVENT_SL,     // speciation and loss
  EVENT_D,      // duplication
  EVENT_T,      // horizontal gene transfer
  EVENT_TL,     // horizontal gene transfer and loss
  EVENT_L,      // loss
  EVENT_None,   // no event
  EVENT_Invalid // invalid event
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

enum class SpeciesTreeAlgorithm {
  User = 0,
  MiniNJ,
  Cherry,
  CherryPro,
  NJst,
  WMinNJ,
  Ustar,
  Random
};

/**
 * Helper methods to work with the enums
 */
class Enums {
public:
  Enums() = delete;

  /**
   * @param m reconciliation model
   * @return the number of free parameters allowed by the model
   */ 
  static unsigned int freeParameters(RecModel m)  {
    switch (m) {
      case RecModel::UndatedDL:
        return 2;
      case RecModel::UndatedDTL:
        return 3;
      case RecModel::UndatedIDTL:
        return 4;
    }
    assert(false);
  }

  /**
   * @param m reconciliation model
   * @return true if the model accounts for horizontal gene transfers
   */
  static bool accountsForTransfers(RecModel m) 
  {
    switch (m) {
    case RecModel::UndatedDL:
      return false;
    case RecModel::UndatedDTL:
    case RecModel::UndatedIDTL:
      return true;
    }
    assert(false);
    return false;
  }

  /**
   * @param m reconciliation model
   * @return true if the corresponding likelihood evaluation
   *         implementation implements a faster approximative
   *         mode (useful to implement heuristics in the search)
   */
  static bool implementsApproxLikelihood(RecModel m)
  {
    switch (m) {
      case RecModel::UndatedDL:
        return false;
      case RecModel::UndatedDTL:
      case RecModel::UndatedIDTL:
        return false;
    }
    assert(false);
    return false;
  }

  static SpeciesTreeAlgorithm strToSpeciesTree(const std::string &str) 
  {
    if (str == std::string("MiniNJ")) {
      return SpeciesTreeAlgorithm::MiniNJ;
    } else if (str == std::string("NJst")) {
      return SpeciesTreeAlgorithm::NJst;
    } else if (str == std::string("WMiniNJ")) {
      return SpeciesTreeAlgorithm::WMinNJ;
    } else if (str == std::string("Ustar")) {
      return SpeciesTreeAlgorithm::Ustar;
    } else if (str == std::string("Cherry")) {
      return SpeciesTreeAlgorithm::Cherry;
    } else if (str == std::string("CherryPro")) {
      return SpeciesTreeAlgorithm::CherryPro;
    } else if (str == std::string("Random") || 
        str == std::string("random")) {
      return SpeciesTreeAlgorithm::Random;
    } else {
      return SpeciesTreeAlgorithm::User;
    }
  }
};


