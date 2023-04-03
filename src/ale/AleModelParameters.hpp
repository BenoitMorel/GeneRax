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
   *  @param cat index to cat label
   *  @param info Description of the model
   */
  AleModelParameters(const Parameters &startingRates,
      const std::vector<unsigned int> &speciesToCategory,
      const std::vector<std::string> &catToLabel,
      const RecModelInfo &info):
    _info(info),
    _speciesToCat(speciesToCategory),
    _catToLabel(catToLabel),
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
  double getRateFromSpecies(unsigned int species, unsigned int rate) const {
    return getRateFromCat(_speciesToCat[species], rate);
  }
  
  double getRateFromCat(unsigned int cat, unsigned int rate) const {
    return _parameters[cat * perCategoryFreeParameters() + rate];
  }

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

  friend std::ostream &operator<<(std::ostream &os, const AleModelParameters &mp)  {
    os << "[";
    for (unsigned int c = 0; c < mp.categoryNumber(); ++c) {
      os << mp._catToLabel[c] << " ";
      os << "(";
      for (unsigned int r = 0; r < mp.perCategoryFreeParameters(); ++r) {
        os << mp.getRateFromCat(c, r);
        if (r !=  mp.perCategoryFreeParameters() -1) {
          os << ",";
        }
      }
      os << ")";
      if (c != mp.categoryNumber() - 1) {
        os << ",";
      }
    }
    os << "]";
    return os;
  }
private:
  RecModelInfo _info;
  std::vector<unsigned int> _speciesToCat;
  std::vector<std::string> _catToLabel;
  unsigned int _catNumber;
  Parameters _parameters;
};



