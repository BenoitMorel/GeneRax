#pragma once

#include <string>
#include <unordered_set>
#include <vector>
#include <util/types.hpp>

class SpeciesTree;
class SpeciesTreeLikelihoodEvaluatorInterface;
class AverageStream;
class SpeciesSearchState;
class Scenario;
static unsigned int hashints(unsigned int a, unsigned int b) {
  return (a + b) * (a + b + 1) / 2  + b;
}

struct TransferMove {
  unsigned int prune;
  unsigned int regraft;
  double transfers;
  TransferMove(): prune(0), regraft(0), transfers(0.0) {}
  TransferMove(unsigned int p, unsigned int r, double t): prune(p), regraft(r), transfers(t) {}
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

namespace std { template<> struct hash<TransferMove>
{
  size_t operator()(const TransferMove & obj) const {
    return hash<int>()(static_cast<int>(
       hashints(hashints(obj.prune, obj.regraft), obj.transfers)));
    }
};}


struct MovesBlackList {
  std::unordered_set<TransferMove> _blacklist;
  bool isBlackListed(const TransferMove &move) { return _blacklist.find(move) != _blacklist.end();}
  void blacklist(const TransferMove &move) { _blacklist.insert(move); }
};

struct PerCorePotentialTransfers {
  
  PerCorePotentialTransfers() {}  
  void addScenario(const Scenario &scenario);
  unsigned int getPotentialTransfers(unsigned int src, unsigned int dest);
  MatrixUint copies; // copies[species][family]
};


class SpeciesTransferSearch {
public:
  
  static bool transferSearch(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    SpeciesSearchState &searchState);

  static void getSortedTransferList(SpeciesTree &speciesTree,
    SpeciesTreeLikelihoodEvaluatorInterface &evaluation,
    unsigned int minTransfers,
    MovesBlackList &blacklist,
    std::vector<TransferMove> &transferMoves
    );

};
