#pragma once

#include <cassert>
#include <string>
#include <unordered_map>
#include <vector>
#include <IO/Logger.hpp>

using StringToUintMap = std::unordered_map<std::string, unsigned int>;

/**
 *  Reconciliation models 
 */
enum class RecModel {
  UndatedDL, 
  UndatedDTL, 
  ParsimonyD,
  SimpleDS
};


/*
 *  DTLRates numerical optimization methods
 */
enum class RecOpt {
  Grid, Simplex, Gradient, None
};


/*
 * Gene tree search mode
 */
enum class GeneSearchStrategy {
  SPR, EVAL, SKIP
};

/**
 * Species tree search mode
 */
enum class SpeciesSearchStrategy {
  SPR, TRANSFERS, HYBRID, REROOT, EVAL, SKIP
};

/**
 *  Transfer constraint
 */
enum class TransferConstaint {
  NONE, PARENTS, SOFTDATED
};

enum class OriginationStrategy {
  UNIFORM, ROOT, LCA
};

/*
 *  Output formats for reconciled gene trees
 */
enum class ReconciliationFormat {
  NHX = 0, RecPhyloXML, NewickEvents
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
 *  Defines how to estimate the root frequencies when
 *  building the conditional clade probabilities from
 *  a list of trees
 */
enum class CCPRooting {
  // the input trees are considered unrooted and all root position
  // have the same frequency
  UNIFORM = 0, 
  // the input trees must be rooted and the alternative 
  // root positions have a null frequency
  ROOTED,
  // the input trees are considered unrooted and MAD rooting
  // is used to estimate the frequency of the root positions
  MAD
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
      case RecModel::ParsimonyD:
        return 0;
      case RecModel::SimpleDS:
        return 1;
      default:
        assert(false);
        return 0;
    }
  }

  static std::vector<std::string> parameterNames(RecModel m)  {
    std::vector<std::string> res;
    switch (m) {
      case RecModel::UndatedDL:
        res.push_back("D");
        res.push_back("L");
        break;
      case RecModel::UndatedDTL:
        res.push_back("D");
        res.push_back("L");
        res.push_back("T");
        break;
      case RecModel::ParsimonyD:
        break;
      case RecModel::SimpleDS:
        res.push_back("D");
        break;
    }
    return res;
  }

  /**
   * @param m reconciliation model
   * @return true if the model accounts for horizontal gene transfers
   */
  static bool accountsForTransfers(RecModel m) 
  {
    switch (m) {
    case RecModel::UndatedDL:
    case RecModel::ParsimonyD:
    case RecModel::SimpleDS:
      return false;
    case RecModel::UndatedDTL:
      return true;
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


  static const char *getEventName(ReconciliationEventType type) {
    static const char *eventNames[]  = {"S", "SL", "D", "T", "TL", "L", "Leaf", "Invalid"};
    return eventNames[static_cast<int>(type)];
  }
  
  static std::string originationToStr(OriginationStrategy os)
  {
    switch(os) {
    case OriginationStrategy::UNIFORM:
      return "UNIFORM";
    case OriginationStrategy::ROOT:
      return "ROOT";
    case OriginationStrategy::LCA:
      return "LCA";
    }
    exit(41);
  }

  static OriginationStrategy strToOrigination(const std::string &str) 
  {
    if (str == "UNIFORM") {
      return OriginationStrategy::UNIFORM;
    }  else if (str == "ROOT") {
      return OriginationStrategy::ROOT;
    }  else if (str == "LCA") {
      return OriginationStrategy::LCA;
    } else {
      Logger::info << "Invalid origination strategy " << str << std::endl;
      exit(41);
    }
  }

};


