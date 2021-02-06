#include <IO/GeneSpeciesMapping.hpp>
#include <trees/PLLRootedTree.hpp>
#include <ccp/ConditionalClades.hpp>

template <class REAL>
class UndatedDLMultiModel {
public: 
  UndatedDLMultiModel(PLLRootedTree &speciesTree, 
      const GeneSpeciesMapping &geneSpeciesMapping, 
      const ConditionalClades &geneDistribution);

  virtual ~UndatedDLMultiModel() {}

  virtual double computeLogLikelihood();

private:
  
  const ConditionalClades &_ccp;
  const PLLRootedTree &_speciesTree; 
  double _PS;
  double _PD;
  double _PL;
  std::vector<double> _uE; // Extinction probability, per species branch
  using DLCLV = std::vector<REAL>;
  std::vector<DLCLV> _dlclvs;


};

template <class REAL>
UndatedDLMultiModel<REAL>::UndatedDLMultiModel(PLLRootedTree &speciesTree, 
    const GeneSpeciesMapping &geneSpeciesMapping, 
    const ConditionalClades &ccp):
  _ccp(ccp),
  _speciesTree(speciesTree),
  _PS(0.7),
  _PD(0.15),
  _PL(0.15),
  _uE(speciesTree.getNodesNumber(), REAL())
{
  std::vector<REAL> zeros(speciesTree.getNodesNumber());
  _dlclvs = std::vector<std::vector<REAL> >(
      _ccp.getCladesNumber(), zeros);

}



template <class REAL>
double UndatedDLMultiModel<REAL>::computeLogLikelihood()
{
  return 0.0;
}

