#pragma once

#include <cassert>

enum class RecModel {
  UndatedDL, UndatedDTL, UndatedDTLAdvanced
};

enum class RecOpt {
  Grid, Simplex, Gradient
};

enum class Strategy {
  SPR, EVAL
};

enum class ReconciliationFormat {
  NHX = 0, RecPhyloXML
};


enum class ReconciliationEventType {
  EVENT_S = 0 , EVENT_SL, EVENT_D, EVENT_T, EVENT_TL, EVENT_L, EVENT_None, EVENT_Invalid
};

class Enums {
public:
  
  static int freeParameters(RecModel m, unsigned int speciesCount = 1)  {
    switch (m) {
      case RecModel::UndatedDL:
        return 2;
      case RecModel::UndatedDTL:
        return 3;
      case RecModel::UndatedDTLAdvanced:
        return 2 + speciesCount;
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
    case RecModel::UndatedDTLAdvanced:
      return true;
    }
    assert(false);
    return false;
  }
};


