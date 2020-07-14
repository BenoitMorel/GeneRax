#pragma once

#include <maths/ModelParameters.hpp>
#include <util/enums.hpp>
#include <util/RecModelInfo.hpp>

/**
 *  Hold the reconciliation rates, information about the model,
 *  and the number of families (relevant is perFamilyRates is set):
 */
class ModelParameters {
public:
  ModelParameters()
  {
  }
  
  /**
   *  @param rates Starting rates 
   *  @param info Reconciliation model information 
   */
  ModelParameters(const Parameters &rates,
      unsigned int familiesNumber,
      const RecModelInfo &info):
    rates((info.perFamilyRates ? familiesNumber : 1), rates),
    info(info),
    familiesNumber(familiesNumber)
  {
    std::cerr << familiesNumber << std::endl;
  
  }

  Parameters getRates(unsigned int familyIndex) const {
    if (info.perFamilyRates) {
      return rates.getSubParameters(familyIndex * info.modelFreeParameters(),
          info.modelFreeParameters());
    } else {
      return rates;
    }
  }
 
  void setRates(unsigned int familyIndex, const Parameters &newRates)
  {
    for (unsigned int i = 0; i < info.modelFreeParameters(); ++i) {
      rates[info.modelFreeParameters() * familyIndex + i] = newRates[i];
    }
  }

public:
  Parameters rates;
  RecModelInfo info;
  unsigned int familiesNumber;
};
