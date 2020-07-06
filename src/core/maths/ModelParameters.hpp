#pragma once

#include <maths/ModelParameters.hpp>
#include <util/enums.hpp>

/**
 *  Hold the reconciliation rates and some information about them:
 *  - The rate values
 *  - Do the families have different rates?
 *  - The reconciliation model (UndatedDL, UndatedDTL etc.)
 *
 *  The number of rates stored is: freeParameters * numberOfRateSets
 *  where freeParameters is the number of free parameters in the 
 *  reconciliation model (3 for UndatedDTL) and numberOfRateSets is
 *  1 if perFamilyRates is false, and familyNumber else.
 */
class ModelParameters {
public:
  ModelParameters():
    model(RecModel::UndatedDL),
    perFamilyRates(false),
    modelFreeParameters(0)
  {
  }
  
  /**
   *  @param rates Starting rates 
   *  @param model Reconciliation model 
   *  @param perFamilyRates Are the rates per family?
   *  @param familyNumber Total number of families
   */
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
