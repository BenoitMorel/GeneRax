#include "ReconciliationBLEstimator.hpp" 

#include <memory>
#include <algorithm>
#include <routines/Routines.hpp>
#include <maths/ModelParameters.hpp>
#include <parallelization/PerCoreGeneTrees.hpp>

static void getAverageDepthRec(pll_rnode_t *node,
    double currentDepth,
    double &sumDepths,
    unsigned int &count)
{
  currentDepth += node->length;
  if (!node->left) {
    count++;
    sumDepths += currentDepth;
  } else {
    getAverageDepthRec(node->left, currentDepth, sumDepths, count);
    getAverageDepthRec(node->right, currentDepth, sumDepths, count);
  }
}

static double getAverageDepth(pll_rnode_t *node)
{
  double sumDepths = 0.0;
  unsigned int count = 0;
  // we do not want to count this node's length
  getAverageDepthRec(node, -node->length, sumDepths, count); 
  assert(count != 0);
  return sumDepths / (double(count));
}

static void balanceRoot(PLLRootedTree &tree) 
{
  auto root = tree.getRoot();
  auto left = root->left;
  auto right = root->right;
  auto initialLength = left->length + right->length;
  auto leftDepth = getAverageDepth(left);
  auto rightDepth = getAverageDepth(right);
  auto diff = leftDepth - rightDepth;
  left->length -= diff / 2.0;
  right->length += diff / 2.0;
  double epsilon = 0.0000001;
  left->length = std::min(std::max(epsilon, left->length), 
      initialLength - epsilon);
  right->length = std::min(std::max(epsilon, right->length), 
      initialLength - epsilon);
}


/**
 *  Fills branch lengths between speciation or leaf events that are direct
 *  parent in the species tree.
 *  Speciation-loss events are not used because they happen along the gene
 *  branch, and not on a gene node, and thus do not hold the time
 *  of speciation
 */
static void estimateBLRecursive(pll_unode_t *node,
    bool isVirtualRoot,
    unsigned int ancestralSpeciesId,
    double lengthToAncestralSpecies,
    pll_rtree_t *speciesTree,
    const std::vector<std::vector<Scenario::Event> > &geneToEvents,
    double familyWeight,
    std::vector<double> &speciesSumBL,
    std::vector<double> &speciesWeightBL)
{
  auto &lastEvent = geneToEvents[node->node_index].back();
  bool isSpeciation = lastEvent.type == ReconciliationEventType::EVENT_S;
  isSpeciation |= lastEvent.type == ReconciliationEventType::EVENT_None;
  lengthToAncestralSpecies += isVirtualRoot ? (node->length / 2.0)
    : node->length;
  if (isSpeciation) {
    auto currentSpeciesId = lastEvent.speciesNode;
    auto currentSpecies = speciesTree->nodes[currentSpeciesId];
    bool isDirectSpeciation = false;
    if (currentSpecies->parent) {
      isDirectSpeciation |= currentSpecies->parent->node_index 
      == ancestralSpeciesId;
    }
    if (isDirectSpeciation) {
      speciesSumBL[currentSpeciesId] += lengthToAncestralSpecies * familyWeight;
      speciesWeightBL[currentSpeciesId] += familyWeight; 
    }
    lengthToAncestralSpecies = 0.0;  
    ancestralSpeciesId = currentSpeciesId;
  }
  if (node->next) {
    // now go down
    auto leftGeneNode = !isVirtualRoot ? node->next->back
      : node->next;
    auto rightGeneNode = !isVirtualRoot ? node->next->next->back
      : node->next->back;
    estimateBLRecursive(leftGeneNode,
        false,
        ancestralSpeciesId,
        lengthToAncestralSpecies,
        speciesTree,
        geneToEvents,
        familyWeight,
        speciesSumBL,
        speciesWeightBL);
    estimateBLRecursive(rightGeneNode,
        false,
        ancestralSpeciesId,
        lengthToAncestralSpecies,
        speciesTree,
        geneToEvents,
        familyWeight,
        speciesSumBL,
        speciesWeightBL);
  }
}

    

void estimateBLForFamily(const Scenario &scenario,
    double familyWeight, 
    std::vector<double> &speciesSumBL,
    std::vector<double> &speciesWeightBL)
{
  auto geneRoot = scenario.getGeneRoot();
  pll_unode_t virtualRoot;
  virtualRoot.next = geneRoot;
  virtualRoot.node_index = scenario.getVirtualRootIndex();
  virtualRoot.label = nullptr;
  virtualRoot.length = 0.0;
  auto ancestralSpeciesId = static_cast<unsigned int>(-1);
  double lengthToAncestralSpecies = 0.0;
  estimateBLRecursive(&virtualRoot,
      true,
      ancestralSpeciesId,
      lengthToAncestralSpecies,
      scenario.getSpeciesTree(),
      scenario.getGeneIdToEvents(),
      familyWeight,
      speciesSumBL,
      speciesWeightBL);
}

double getFamilyWeight(const FamilyInfo &info) 
{
  if (info.alignmentFile == "" || info.libpllModel == "") {
    return 1.0;
  }
  auto alignmentLength = LibpllParsers::getMSAEntropy(info.alignmentFile,
      info.libpllModel);
  if (alignmentLength == 0) {
    return 1.0;
  }
  return static_cast<double>(alignmentLength);
}

void ReconciliationBLEstimator::estimate(
      const std::string &speciesTreeFile,
      const Families &families, 
      const RecModelInfo &recModelInfo,
      const Parameters &rates)
{
    PLLRootedTree speciesTree(speciesTreeFile, true);
    PerCoreGeneTrees geneTrees(families);
    ModelParameters modelParameters(rates, 
        families.size(),
        recModelInfo);
    const unsigned int samples = 0;
    const bool optimizeRates = false;
    std::vector<Scenario> scenarios;
    Logger::timed << "[BL estimation] Infering scenarios" << std::endl;
    Routines::inferAndGetReconciliationScenarios(speciesTree,
        geneTrees,
        modelParameters,
        samples,
        optimizeRates,
        scenarios);
    double speciesNodesNumber = speciesTree.getNodesNumber();
    std::vector<double> speciesSumBL(speciesNodesNumber, 0.0);
    std::vector<double> speciesWeightBL(speciesNodesNumber, 0.0);
       
    for (unsigned int i = 0; i < geneTrees.getTrees().size(); ++i) {
      double familyWeight = getFamilyWeight(families[geneTrees.getTrees()[i].familyIndex]);
      estimateBLForFamily(scenarios[i],
        familyWeight,
        speciesSumBL,
        speciesWeightBL);
    }
    ParallelContext::sumVectorDouble(speciesSumBL);
    ParallelContext::sumVectorDouble(speciesWeightBL);
    for (unsigned int i = 0; i < speciesSumBL.size(); ++i) {
      auto length = 0.0;
      if (speciesWeightBL[i] != 0.0) {
        length = speciesSumBL[i] / speciesWeightBL[i];
      }
      speciesTree.getNode(i)->length = length;
    }
    balanceRoot(speciesTree);
    Logger::timed << "[BL estimation] Inferred branch lengths:\n" << speciesTree.getNewickString() << std::endl;
    if (ParallelContext::getRank() == 0) {
      speciesTree.save(speciesTreeFile);
    }
    Logger::timed << "[BL estimation] Done" << std::endl;
}

