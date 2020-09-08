#include <string>
#include <IO/Families.hpp>
class RecModelInfo;
class Parameters;



class ReconciliationBLEstimator {
public:
  ReconciliationBLEstimator() = delete;

  static void estimate(const std::string &speciesTreeFile,
      const Families &families, 
      const RecModelInfo &recModelInfo,
      const Parameters &rates);

  




};


