#pragma once

#include <cassert>

enum RecModel {
  UndatedDL, UndatedDTL
};

class Enums {
public:
  
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
  Grid, Simplex
};

enum Strategy {
  SPR, EVAL
};

