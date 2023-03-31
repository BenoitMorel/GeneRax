#pragma once

#include <util/RecModelInfo.hpp>

#include <vector>
#include <algorithm>

class AleModelParameters {
public:
  AleModelParameters():
    _catNumber(0)
  {

  }
  
  AleModelParameters(const Parameters &startingRates,
      const std::vector<unsigned int> &speciesToCategory,
      const RecModelInfo &info):
    _info(info),
    _speciesToCat(speciesToCategory),
    _catNumber(*std::max_element(std::begin(speciesToCategory), std::end(speciesToCategory)) + 1),
    _parameters(_catNumber, startingRates)
  {
    for (unsigned int i = 0; i < getTotalFreeParameters(); ++i) {
      _parameters[i] = startingRates[i % perCategoryFreeParameters()];
    }
  }

  unsigned int categoryNumber() const {return _catNumber;}

  unsigned int perCategoryFreeParameters() const {return _info.modelFreeParameters();}
  
  unsigned int getTotalFreeParameters() const {return categoryNumber() * perCategoryFreeParameters();} 
  
  double getRate(unsigned int species, unsigned int rate) const {
    return _parameters[_speciesToCat[species] * perCategoryFreeParameters() + rate];}

  const RecModelInfo &getInfo() const {return _info;}

  void setRates(const Parameters &parameters) {
    assert(parameters.dimensions() == getTotalFreeParameters());
    _parameters = parameters;
  }

  const Parameters &getRates() {return _parameters;}

private:
  RecModelInfo _info;
  std::vector<unsigned int> _speciesToCat;
  unsigned int _catNumber;
  Parameters _parameters;
};



