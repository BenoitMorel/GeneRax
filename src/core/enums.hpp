#pragma once

#include <cassert>

enum RecModel {
  UndatedDL, UndatedDTL, DatedDL
};

class Enums {
public:
  
  static bool accountsForTransfers(RecModel m) 
  {
    switch (m) {
    case UndatedDL: case DatedDL:
      return false;
    case UndatedDTL:
      return true;
    }
    assert(false);
  }
};

enum RecOpt {
  Grid, Simplex
};

enum Strategy {
  SPR, EVAL
};

