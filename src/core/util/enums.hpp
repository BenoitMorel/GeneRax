#pragma once

#include <cassert>

enum RecModel {
  UndatedDL, UndatedDTL
};

class Enums {
public:
  
  static int freeParameters(RecModel m)  {
    return m == UndatedDL ? 2 : 3;
  }

  static bool accountsForTransfers(RecModel m) 
  {
    switch (m) {
    case UndatedDL:
      return false;
    case UndatedDTL:
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
  SIMPLE_SEARCH, SUBSAMPLE_SEARCH
};

enum ReconciliationFormat {
  NHX = 0, RecPhyloXML
};


enum ReconciliationEventType {
  EVENT_S = 0 , EVENT_SL, EVENT_D, EVENT_T, EVENT_TL, EVENT_L, EVENT_None, EVENT_Invalid
};

