#include "Families.hpp"
#include <likelihoods/LibpllEvaluation.hpp>
#include <IO/Logger.hpp>

bool filterFamily(const FamilyInfo &family)
{
  bool res = true;
  try {
    res &= LibpllEvaluation::isValidAlignment(family.alignmentFile, family.libpllModel);
  } catch (...) {
    res = false;
  }
  return res;
}

void filterFamilies(Families &families)
{
  Families copy = families;
  families.clear();
  unsigned int invalid = 0;
  for (auto &family: copy) {
    if (filterFamily(family)) {
      families.push_back(family);
    } else {
      invalid++;
    }
  }
  if (invalid) {
    Logger::error << "Found " << invalid << " invalid families (they will be discarded from the analysis)" << std::endl;
  }   
}



