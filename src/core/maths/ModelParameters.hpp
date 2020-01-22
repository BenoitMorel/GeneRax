#pragma once

#include <maths/ModelParameters.hpp>
#include <util/enums.hpp>

class ModelParameters {
public:


  ModelParameters(const Parameters &rates, 
      RecModel model,
      bool perFamilyRates,
      unsigned int familyNumber):
    rates((perFamilyRates ? familyNumber : 1), rates),
    model(model),
    perFamilyRates(perFamilyRates),
    modelFreeParameters(Enums::freeParameters(model))
  {}

  Parameters getRates(unsigned int familyIndex) const {
    if (perFamilyRates) {
      return rates.getSubParameters(familyIndex * modelFreeParameters, modelFreeParameters);  
    } else {
      return rates;
    }
  }
 
  void setRates(unsigned int familyIndex, const Parameters &newRates)
  {
    for (unsigned int i = 0; i < modelFreeParameters; ++i) {
      rates[modelFreeParameters * familyIndex + i] = newRates[i];
    }
  }

public:
  Parameters rates;
  RecModel model;
  bool perFamilyRates;
private:
  unsigned int modelFreeParameters;
};
