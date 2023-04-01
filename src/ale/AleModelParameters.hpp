#pragma once

#include <util/RecModelInfo.hpp>

#include <vector>
#include <algorithm>
#include <fstream>

/**
 * Stores the model parameters: DTL rates, global or per species category
 * A category is a set of species that shares the same DTL rates
 */
class AleModelParameters {
public:
  /**
   *  @brief Default constructor
   */
  AleModelParameters():
    _catNumber(0)
  {
    
  }
 
  /**
   *  @brief Constructor
   *  @param startingRates the set of rates to assign to each category
   *  @param speciesToCategory mapping species node index to category
   *  @param info Description of the model
   */
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
  
  /**
   *  Number of species categories
   */
  unsigned int categoryNumber() const {return _catNumber;}

  /**
   *  Number of free parameters for a single species category
   */
  unsigned int perCategoryFreeParameters() const {return _info.modelFreeParameters();}
  
  /**
   * Number of tree parameters for all categories
   */
  unsigned int getTotalFreeParameters() const {return categoryNumber() * perCategoryFreeParameters();} 
  
  /**
   *  Get a specific rate for a specific species
   */
  double getRate(unsigned int species, unsigned int rate) const {
    return _parameters[_speciesToCat[species] * perCategoryFreeParameters() + rate];}

  /**
   *  Get model description
   */
  const RecModelInfo &getInfo() const {return _info;}

  /**
   *  Set rates from a parameter vector
   */
  void setRates(const Parameters &parameters) {
    assert(parameters.dimensions() == getTotalFreeParameters());
    _parameters = parameters;
  }

  /**
   *  Get rates as a parameter vector
   */
  const Parameters &getRates() {return _parameters;}

  friend std::ofstream &operator<<(std::ofstream &os, const ModelParameters &mp)  {
    os << "[";
    for (unsigned int c = 0; c < categoryNumber() << ++c) {
      os << "(";
      for (unsigned int r = 0; r < perCategoryFreeParameters(); ++r) {
           
      }
      os << ")";
      if (c != categoryNumber() - 1) {
        os << ",";
      }
    }
    os << "]";
i   return os;
  }
private:
  RecModelInfo _info;
  std::vector<unsigned int> _speciesToCat;
  unsigned int _catNumber;
  Parameters _parameters;
};



