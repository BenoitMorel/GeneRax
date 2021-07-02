#include "SpeciesTransferSearch.hpp"

#include <algorithm>

#include <search/SpeciesSearchCommon.hpp>
#include <trees/SpeciesTree.hpp>
#include <trees/PLLRootedTree.hpp>

struct TransferMove {
  unsigned int prune;
  unsigned int regraft;
  double transfers;
  TransferMove(): prune(0), regraft(0), transfers(0.0) {
  }
  TransferMove(unsigned int p, unsigned int r, double t): prune(p), regraft(r), transfers(t) {
  }
  bool operator < (const TransferMove& tm) const
  {
    if (transfers != tm.transfers) {
      return transfers > tm.transfers;
    } else if (regraft != tm.regraft) {
      return regraft > tm.regraft;
    } else {
      return prune > tm.prune;
    }
  }
  bool operator ==(const TransferMove& obj) const
  {
    return (obj.prune == prune) && (obj.regraft == regraft) && (obj.transfers == transfers); 
  }
};

static unsigned int hashints(unsigned int a, unsigned int b)
{
  return (a + b) * (a + b + 1) / 2  + b;
}

namespace std {
template<>
struct hash<TransferMove>
{
  size_t
    operator()(const TransferMove & obj) const
    {
      return hash<int>()(static_cast<int>(
            hashints(hashints(obj.prune, obj.regraft), obj.transfers)));
    }
};
}

struct MovesBlackList {
  std::unordered_set<TransferMove> _blacklist;
  bool isBlackListed(const TransferMove &move) { return _blacklist.find(move) != _blacklist.end();}
  void blacklist(const TransferMove &move) { _blacklist.insert(move); }
};


static bool transferRound(SpeciesTree &speciesTree,
  SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
  SpeciesSearchState &searchState,
  double previousBestLL,
  double &newBestLL,
  MovesBlackList &blacklist,
  bool &maxImprovementsReached)
{
  newBestLL = previousBestLL;
  maxImprovementsReached = false;
  unsigned int minTransfers = 1;
  auto hash1 = speciesTree.getNodeIndexHash(); 
  TransferFrequencies frequencies;
  PerSpeciesEvents perSpeciesEvents;
  evaluation.getTransferInformation(speciesTree.getTree(),
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
  std::vector<TransferMove> transferMoves;
  for (unsigned int from = 0; from < frequencies.count.size(); ++from) {
    for (unsigned int to = 0; to < frequencies.count.size(); ++to) {
      auto regraft = labelsToIds[frequencies.idToLabel[from]];
      auto prune = labelsToIds[frequencies.idToLabel[to]];
      auto count = frequencies.count[from][to];
      transfers += count;
      if (count < minTransfers) {
        continue;
      }
      if (SpeciesTreeOperator::canApplySPRMove(speciesTree, prune, regraft)) {
        TransferMove move(prune, regraft, count);
        double factor = 1.0;
        if (evaluation.pruneSpeciesTree()) {
          factor /= (1.0 + sqrt(speciesFrequencies[prune]));
          factor /= (1.0 + sqrt(speciesFrequencies[regraft]));
        }
        if (!blacklist.isBlackListed(move)) {
          transferMoves.push_back(TransferMove(prune, regraft, factor * count)); 
        }
      }
    }
  }
  Logger::timed << "[Species search] Start new transfer-guided round. Inferred transfers:" << transfers << std::endl;
  std::sort(transferMoves.begin(), transferMoves.end());
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
      double testedLL = 0.0;
      if (SpeciesSearchCommon::testSPR(speciesTree, evaluation, 
            transferMove.prune, transferMove.regraft, searchState, 
            newBestLL, testedLL)) {
        newBestLL = testedLL;
        failures = 0;
        improvements++;
        alreadyPruned.insert(transferMove.prune);
        auto pruneNode = speciesTree.getNode(transferMove.prune);
        Logger::timed << "\tbetter tree (transfers:" 
          << transferMove.transfers 
          << ", trial: " 
          << trials 
          << ", ll=" 
          << newBestLL 
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
          if (SpeciesSearchCommon::veryLocalSearch(speciesTree,
              evaluation,
              searchState,
              transferMove.prune,
              newBestLL,
              testedLL)) {
            newBestLL = testedLL;
          }
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
        return improvements > 1;
      }
    }  
  }
  return improvements > 1;
}


bool SpeciesTransferSearch::transferSearch(
  SpeciesTree &speciesTree,
  SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
  SpeciesSearchState &searchState,
  double previousBestLL,
  double &newBestLL)
{
  newBestLL = previousBestLL;
  Logger::info << std::endl;
  Logger::timed << "[Species search]" 
    << " Starting species tree transfer-guided search, bestLL=" << previousBestLL 
    <<  ", hash=" << speciesTree.getHash() << ")" << std::endl;
  MovesBlackList blacklist;
  bool maxImprovementsReached = true;
  bool stop = false;
  bool better = false;
  while (!stop) {
    if (maxImprovementsReached) {
      Logger::info << "TODO: OPTIMIZE RATES" << std::endl;
      //newBestLL = hack.optimizeDTLRates();
    }
    stop = !transferRound(speciesTree, evaluation, searchState, 
        previousBestLL, newBestLL, blacklist, maxImprovementsReached);
    if (!stop) {
      better = true;
    }
  } 
  Logger::timed << "[Species search] After transfer search: LL=" << newBestLL << std::endl;
  return better;
}


