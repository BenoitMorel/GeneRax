#pragma once

#include <cassert>

enum class RecModel {
  UndatedDL, UndatedDTL
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
};


