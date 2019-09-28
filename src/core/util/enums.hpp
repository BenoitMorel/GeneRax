#pragma once

#include <cassert>

enum RecModel {
  UndatedDL, UndatedDTL, UndatedDTLAdvanced
};

class Enums {
public:
  
  static int freeParameters(RecModel m, unsigned int speciesCount = 1)  {
    switch (m) {
      case UndatedDL:
        return 2;
      case UndatedDTL:
        return 3;
      case UndatedDTLAdvanced:
        return 2 + speciesCount;
    }
    assert(false);
  }

  static bool accountsForTransfers(RecModel m) 
  {
    switch (m) {
    case UndatedDL:
      return false;
    case UndatedDTL:
      return true;
    case UndatedDTLAdvanced:
      return true;
    }
    assert(false);
    return false;
  }
};

enum RecOpt {
  Grid, Simplex, Gradient
};

enum Strategy {
  SPR, EVAL
};

enum SpeciesRaxStrategy {
  SIMPLE_SEARCH
};

enum ReconciliationFormat {
  NHX = 0, RecPhyloXML
};


enum ReconciliationEventType {
  EVENT_S = 0 , EVENT_SL, EVENT_D, EVENT_T, EVENT_TL, EVENT_L, EVENT_None, EVENT_Invalid
};

