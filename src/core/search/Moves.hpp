#pragma once

#include <likelihoods/LibpllEvaluation.hpp>
#include <search/Rollbacks.hpp>

#include <memory>



class JointTree;


class SPRMove {
public:
  SPRMove(unsigned int pruneIndex, unsigned int regraftIndex, const std::vector<unsigned int> &path);
  virtual ~SPRMove() {}
  static std::shared_ptr<SPRMove> createSPRMove(unsigned int pruneIndex, 
      unsigned int regraftIndex, 
      const std::vector<unsigned int> &path);
  std::shared_ptr<SPRRollback> applyMove(JointTree &tree);
  void optimizeMove(JointTree &tree);
  std::ostream& print(std::ostream & os) const;
  void setScore(double score) {_score = score;}
  double getScore() const {return _score;}
  
private:
  unsigned int pruneIndex_;
  unsigned int regraftIndex_;
  std::vector<unsigned int> path_;
  std::vector<corax_unode_t *> branchesToOptimize_;
  double _score;
};

