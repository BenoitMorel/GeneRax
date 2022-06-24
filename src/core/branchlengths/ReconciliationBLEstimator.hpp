#include <string>
#include <IO/Families.hpp>
class ModelParameters;



class ReconciliationBLEstimator {
public:
  ReconciliationBLEstimator() = delete;

  /**
   *  Estimate the branch lengths of a species tree by 
   *  reconciliating the gene trees (in families) using
   *  the model of reconciliation described in modelParameters
   *  and then averaging over the gene tree branch lengths
   *  that represent a branch in the species tree.
   *  The exact procedure is formally described in the SpeciesRax
   *  paper.
   */
  static void estimate(const std::string &speciesTreeFile,
      const Families &families, 
      const ModelParameters &modelParameters);
};


