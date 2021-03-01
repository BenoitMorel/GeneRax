#include <cstdio> 
#include <cstdlib> 
#include <cassert>
#include <iostream>
#include <treetests/AUTest.hpp>

int main(int, char**)
{
  std::cout << "Testing au test..." << std::endl;
  size_t sites = 1000;
  size_t alignments = 3;
  size_t nboot = 10000; 
  std::vector<Likelihoods> likelihoods(alignments);
  for (unsigned int k = 0; k < alignments; ++k) {
    likelihoods[k].resize(sites);
  }
  for (unsigned int i = 0; i < sites; ++i) {
    likelihoods[0][i] = -1.0; 
    if (i < sites / 2) {
      likelihoods[1][i] = -1.001; 
    } else {
      likelihoods[1][i] = -0.999; 
    }
    likelihoods[2][i] = -1.01;
  }

  std::vector<double> pvalues;
  performAUTest(likelihoods, nboot, pvalues); 
  return 0;
}


