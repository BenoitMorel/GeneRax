#include "SpeciesTransferSearch.hpp"

#include <algorithm>

#include <search/SpeciesSearchCommon.hpp>
#include <trees/SpeciesTree.hpp>
#include <trees/PLLRootedTree.hpp>
#include <IO/FileSystem.hpp>



void SpeciesTransferSearch::getSortedTransferList(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    unsigned int minTransfers,
    MovesBlackList &blacklist,
    std::vector<TransferMove> &transferMoves
    )
{
  TransferFrequencies frequencies;
  PerSpeciesEvents perSpeciesEvents;
  evaluation.getTransferInformation(speciesTree,
      frequencies,
      perSpeciesEvents);
  unsigned int speciesNumber = speciesTree.getTree().getNodesNumber();
  std::vector<double> speciesFrequencies;
  for (unsigned int e = 0; e < speciesNumber; ++e) {
    auto &speciesEvents = perSpeciesEvents.events[e];
    speciesFrequencies.push_back(speciesEvents.speciesFrequency());
  }
  unsigned int transfers = 0;
  ParallelContext::barrier();
  std::unordered_map<std::string, unsigned int> labelsToIds;
  speciesTree.getLabelsToId(labelsToIds);
  for (unsigned int from = 0; from < frequencies.count.size(); ++from) {
    for (unsigned int to = 0; to < frequencies.count.size(); ++to) {
      auto regraft = labelsToIds[frequencies.idToLabel[from]];
      auto prune = labelsToIds[frequencies.idToLabel[to]];
      auto count = frequencies.count[from][to];
      transfers += count;
      if (count < minTransfers) {
        continue;
      }
      TransferMove move(prune, regraft, count);
      double factor = 1.0;
      if (true) {//evaluation.pruneSpeciesTree()) {
        factor /= (1.0 + sqrt(speciesFrequencies[prune]));
        factor /= (1.0 + sqrt(speciesFrequencies[regraft]));
      }
      if (!blacklist.isBlackListed(move)) {
        transferMoves.push_back(TransferMove(prune, regraft, factor * count)); 
      }
    }
  }
  Logger::timed << "[Species search] Start new transfer-guided round. Inferred transfers:" << transfers << std::endl;
  std::sort(transferMoves.begin(), transferMoves.end());

}

static bool transferRound(SpeciesTree &speciesTree,
  SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
  SpeciesSearchState &searchState,
  MovesBlackList &blacklist,
  bool &maxImprovementsReached)
{
  maxImprovementsReached = false;
  auto hash1 = speciesTree.getNodeIndexHash(); 
  unsigned int minTransfers = 1;
  std::vector<TransferMove> transferMoves;
 
  SpeciesTransferSearch::getSortedTransferList(speciesTree, evaluation, minTransfers, blacklist, transferMoves);
  auto copyTransferMoves = transferMoves;
  transferMoves.clear();
  for (const auto &transferMove: copyTransferMoves) {
    if (SpeciesTreeOperator::canApplySPRMove(speciesTree, transferMove.prune, transferMove.regraft)) {
      transferMoves.push_back(transferMove); 
    }
  }
  unsigned int speciesNumber = speciesTree.getTree().getNodesNumber();
  unsigned int index = 0;
  const unsigned int stopAfterFailures = 50u;
  const unsigned int stopAfterImprovements = std::max(15u, speciesNumber / 4);
  const unsigned int minTrial = std::max(50u, speciesNumber / 2);
  unsigned int failures = 0;
  unsigned int improvements = 0;
  std::unordered_set<unsigned int> alreadyPruned;
  unsigned int trials = 0;
  for (auto &transferMove: transferMoves) {
    index++;
    if (alreadyPruned.find(transferMove.prune) != alreadyPruned.end()) {
      continue;
    }
    if (SpeciesTreeOperator::canApplySPRMove(speciesTree, transferMove.prune, transferMove.regraft)) {
      blacklist.blacklist(transferMove);
      trials++;
      /*
      Logger::info << "Test "   
          << speciesTree.getNode(transferMove.prune)->label
          << " -> " 
          << speciesTree.getNode(transferMove.regraft)->label
          << std::endl;
          */
      if (SpeciesSearchCommon::testSPR(speciesTree, evaluation, 
            searchState, transferMove.prune, transferMove.regraft)) {
        failures = 0;
        improvements++;
        alreadyPruned.insert(transferMove.prune);
        auto pruneNode = speciesTree.getNode(transferMove.prune);
        Logger::timed << "\tbetter tree (transfers:" 
          << transferMove.transfers 
          << ", trial: " 
          << trials 
          << ", ll=" 
          << searchState.bestLL 
          << ", hash=" 
          << speciesTree.getHash() 
          << ") "  
          << pruneNode->label 
          << " -> " 
          << speciesTree.getNode(transferMove.regraft)->label
          << std::endl;
        // we enough improvements to recompute the new transfers
        hash1 = speciesTree.getNodeIndexHash(); 
        assert(ParallelContext::isIntEqual(hash1));
        if (!searchState.farFromPlausible) {
          SpeciesSearchCommon::veryLocalSearch(speciesTree,
              evaluation,
              searchState,
              transferMove.prune);
        }
      } else {
        failures++;
      }
      bool stop = index > minTrial && failures > stopAfterFailures;
      maxImprovementsReached = improvements > stopAfterImprovements;
      stop |= maxImprovementsReached;
      if (stop) {
        if (searchState.farFromPlausible && !maxImprovementsReached) {
          Logger::timed << 
            "[Species search] Switch to hardToFindBetter mode" << std::endl;
          searchState.farFromPlausible = false;
        }
        return improvements > 0;
      }
    }  
  }
  return improvements > 0;
}


bool SpeciesTransferSearch::transferSearch(
  SpeciesTree &speciesTree,
  SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
  SpeciesSearchState &searchState)
{
  Logger::info << std::endl;
  Logger::timed << "[Species search]" 
    << " Starting species tree transfer-guided search, bestLL=" << searchState.bestLL
    <<  ", hash=" << speciesTree.getHash() << ")" << std::endl;
  MovesBlackList blacklist;
  bool maxImprovementsReached = true;
  bool stop = false;
  bool better = false;
  while (!stop) {
    searchState.bestLL = evaluation.optimizeModelRates();
    stop = !transferRound(speciesTree, evaluation, searchState, 
        blacklist, maxImprovementsReached);
    if (!stop) {
      better = true;
    }
  } 
  Logger::timed << "[Species search] After transfer search: LL=" << searchState.bestLL << std::endl;
  return better;
}


