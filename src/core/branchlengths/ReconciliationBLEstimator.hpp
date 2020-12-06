#include <string>
#include <IO/Families.hpp>
class ModelParameters;



class ReconciliationBLEstimator {
public:
  ReconciliationBLEstimator() = delete;

  static void estimate(const std::string &speciesTreeFile,
      const Families &families, 
      const ModelParameters &parameters);

  




};


