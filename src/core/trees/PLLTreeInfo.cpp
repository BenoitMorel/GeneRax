#include "PLLTreeInfo.hpp"

#include <sstream>
#include <fstream>
#include <IO/Logger.hpp>
#include <maths/Random.hpp>
#include <IO/LibpllException.hpp>

#include <corax/corax.h>
const double DEFAULT_BL = 0.1;

static unsigned int getBestLibpllAttribute() {
  corax_hardware_probe();
  unsigned int arch = CORAX_ATTRIB_ARCH_CPU;
  if (corax_hardware.avx2_present) {
    arch = CORAX_ATTRIB_ARCH_AVX2;
  } else if (corax_hardware.avx_present) {
    arch = CORAX_ATTRIB_ARCH_AVX;
  } else if (corax_hardware.sse_present) {
    arch = CORAX_ATTRIB_ARCH_SSE;
  }
  arch |= CORAX_ATTRIB_SITE_REPEATS;
  return  arch;
}


void treeinfoDestroy(corax_treeinfo_t *treeinfo)
{
  if (!treeinfo)
    return;
  corax_partition_destroy(treeinfo->partitions[0]);
  corax_treeinfo_destroy(treeinfo);
}


PLLTreeInfo::PLLTreeInfo(const std::string &newickStrOrFile,
    bool isNewickAFile,
    const std::string& alignmentFilename,
    const std::string &modelStrOrFile) :
  _treeinfo(nullptr, treeinfoDestroy),
  _model(LibpllParsers::getModel(modelStrOrFile))
{
  std::cerr << "create tree info from " << newickStrOrFile << std::endl;
  PLLSequencePtrs sequences;
  unsigned int *patternWeights = nullptr;
  LibpllParsers::parseMSA(alignmentFilename, _model->charmap(), sequences, patternWeights);
  buildTree(newickStrOrFile, isNewickAFile, sequences);
  auto partition = buildPartition(sequences, patternWeights);
  _treeinfo = std::unique_ptr<corax_treeinfo_t, void(*)(corax_treeinfo_t*)>(
      buildTreeInfo(*_model, partition, *_utree), treeinfoDestroy);
  free(patternWeights);
}
  


void PLLTreeInfo::buildModel(const std::string &modelStrOrFile)
{
  std::string modelStr = modelStrOrFile;
  std::ifstream f(modelStr);
  if (f.good()) {
    getline(f, modelStr);
    modelStr = modelStr.substr(0, modelStr.find(","));
  }
  _model = std::make_unique<Model>(modelStr);
  assert(_model->num_submodels() == 1);

}
void PLLTreeInfo::buildTree(const std::string &newickStrOrFile, 
      bool isNewickAFile, 
      const PLLSequencePtrs &sequences)
{
  if (newickStrOrFile == "__random__") {
    std::vector<const char*> labels;
    for (const auto &seq: sequences) {
      labels.push_back(seq->label);
    }
    unsigned int seed = static_cast<unsigned int>(Random::getInt());
    _utree = std::unique_ptr<PLLUnrootedTree>(new PLLUnrootedTree(labels, seed));
  } else {
    _utree = std::unique_ptr<PLLUnrootedTree>(new PLLUnrootedTree(newickStrOrFile, isNewickAFile));
  }
  _utree->setMissingBranchLengths(DEFAULT_BL);
}

corax_partition_t * PLLTreeInfo::buildPartition(const PLLSequencePtrs &sequences, 
  unsigned int *patternWeights)  
{  
  unsigned int attribute = getBestLibpllAttribute();
  unsigned int tipNumber = static_cast<unsigned int>(sequences.size());
  unsigned int innerNumber = tipNumber -1;
  unsigned int edgesNumber = 2 * tipNumber - 1;
  unsigned int sitesNumber = sequences[0]->len;
  unsigned int ratesMatrices = 1;
  corax_partition_t *partition = corax_partition_create(tipNumber,
      innerNumber,
      _model->num_states(),
      sitesNumber, 
      ratesMatrices, 
      edgesNumber,// prob_matrices
      _model->num_ratecats(),  
      edgesNumber,// scalers
      attribute);  
  if (!partition) 
    throw LibpllException("Could not create libpll partition");
  corax_set_pattern_weights(partition, patternWeights);
  
  // fill partition
  std::map<std::string, unsigned int> tipsLabelling;
  unsigned int labelIndex = 0;
  for (auto &seq: sequences) {
    tipsLabelling[seq->label] = labelIndex;
    corax_set_tip_states(partition, labelIndex, _model->charmap(), seq->seq);
    labelIndex++;
  }
  assign(partition, *_model);
  
  // map tree to partition
  for (auto leaf: _utree->getLeaves()) {
    assert(tipsLabelling.find(leaf->label) != tipsLabelling.end());
    leaf->clv_index = tipsLabelling[leaf->label];
  }
  return partition;
}

corax_treeinfo_t *PLLTreeInfo::buildTreeInfo(const Model &model,
    corax_partition_t *partition,
    PLLUnrootedTree &utree)
{
  // treeinfo
  int params_to_optimize = model.params_to_optimize();
  params_to_optimize |= CORAX_OPT_PARAM_BRANCHES_ITERATIVE;
  std::vector<unsigned int> params_indices(model.num_ratecats(), 0); 
  auto treeinfo = corax_treeinfo_create(utree.getAnyInnerNode(), 
      utree.getLeafNumber(), 1, CORAX_BRLEN_SCALED);
  if (!treeinfo || !treeinfo->root)
    throw LibpllException("Cannot create treeinfo");
  corax_treeinfo_init_partition(treeinfo, 0, partition,
      params_to_optimize,
      model.gamma_mode(),
      model.alpha(), 
      &params_indices[0],
      model.submodel(0).rate_sym().data());
  return treeinfo;
}
