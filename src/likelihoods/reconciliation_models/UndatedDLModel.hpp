#ifndef JOINTSEARCH_UNDATEDDLMODEL_HPP_
#define JOINTSEARCH_UNDATEDDLMODEL_HPP_

#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>
#include <likelihoods/LibpllEvaluation.hpp>
#include <parsers/GeneSpeciesMapping.hpp>
#include <iostream>
#include <limits>
//#define JS_SCALE_FACTOR 115792089000000000000.0  /*  2**256 (exactly)  */
#define JS_SCALE_FACTOR 115792089237316195423570985008687907853269984665640564039457584007913129639936.0  /*  2**256 (exactly)  */
#define JS_SCALE_THRESHOLD (1.0/JS_SCALE_FACTOR)

using namespace std;

#define BIG 100000000

struct ScaledValue {
  ScaledValue():value(0), scaler(BIG) {

  }
  ScaledValue(double v, int s):value(v), scaler(s) {
    if (value == 0) {
      scaler = BIG;
      value = 0;
    }
    else {
      while (value < JS_SCALE_THRESHOLD) {
        scaler += 1;
        value *= JS_SCALE_FACTOR;
      }
    }
  } 

  void add(ScaledValue v) {
    if (v.scaler == scaler) {
      value += v.value;
    } else if (v.scaler != BIG && v.scaler < scaler) {
      scaler = v.scaler;
      value = v.value;
    }
  }

  double getLogValue() {
    if (scaler == BIG) {
      return -std::numeric_limits<double>::infinity();
    }
    return log(value) + scaler * log(JS_SCALE_THRESHOLD);
  }

  double value;
  int scaler;
};

struct ListScaledValue {
  void add(double value, int scaler) {
    add(ScaledValue(value, scaler));
  }
  void add(ScaledValue value) {
    total.add(value);     
  }
  ScaledValue total;
};

class UndatedDLModel: public AbstractReconciliationModel {
public:

  // model
  vector<double> PD; // Duplication probability, per branch
  vector<double> PL; // Loss probability, per branch
  vector<double> PS; // Speciation probability, per branch
  const double O_R; // what is this?

  // SPECIES
  vector<double> uE; // Probability for a gene to become extinct on each brance
  
  // CLVs
  vector<vector<double> > uq;
  vector<vector<int> > uq_scalers;
  vector<double> ll;
  vector<int> ll_scalers;
  
  vector<int> geneIds;
public:
  UndatedDLModel();
  virtual ~UndatedDLModel();
  
  // unherited from parents
  virtual void setRates(double dupRate, double lossRate, double transferRate = 0.0);
  virtual double computeLogLikelihood(shared_ptr<pllmod_treeinfo_t> treeinfo);
  
private:
  void computeProbability(pll_unode_t *geneNode, pll_rnode_t *speciesNode, 
      ScaledValue &proba,
      bool isVirtualRoot = false);
  void updateCLV(pll_unode_t *geneNode);
  void updateCLVs(pllmod_treeinfo_t &treeinfo);
  pll_unode_t *computeLikelihoods(pllmod_treeinfo_t &treeinfo);

};

#endif

