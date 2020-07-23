#pragma once

#include <util/enums.hpp>

struct RecModelInfo {
  RecModel model;
  bool perFamilyRates;

  RecModelInfo():
    model(RecModel::UndatedDTL),
    perFamilyRates(false)
  {

  }

  RecModelInfo(RecModel model,
      bool perFamilyRates):
    model(model),
    perFamilyRates(perFamilyRates)
  {

  }

  unsigned int modelFreeParameters() const {
    return Enums::freeParameters(model);
  }

};
