#pragma once


class SpeciesSplitScore {
public:
  virtual ~SpeciesSplitScore() {};
  virtual void updateSpeciesTree(PLLRootedTree &speciesTree) = 0;
  virtual double getScore() = 0;
};
  

