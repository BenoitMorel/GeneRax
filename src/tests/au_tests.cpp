#include <cstdio> 
#include <cstdlib> 
#include <cassert>
#include <iostream>
#include <treetests/AUTest.hpp>

double randDouble(double max)
{
  return (max / double(RAND_MAX)) * double(rand());
}

int main(int, char**)
{
  srand(42);
  std::cout << "Testing au test..." << std::endl;
  size_t sites = 1000;
  size_t alignments = 3;
  size_t nboot = 10000; 
  std::vector<Likelihoods> likelihoods(alignments);
  for (unsigned int k = 0; k < alignments; ++k) {
    likelihoods[k].resize(sites);
  }
  for (unsigned int i = 0; i < sites; ++i) {
    likelihoods[0][i] = -randDouble(1.0); 
    likelihoods[1][i] = -randDouble(1.1); //(double(rand() % 1500));
    likelihoods[2][i] = -0.5; //(double(rand() % 1500));
    for (unsigned int k = 0; k < alignments; k++) {
      likelihoods[k][i] *= 1000.0;
    }
  }

  std::vector<double> pvalues;
  performAUTest(likelihoods, nboot, pvalues); 
  return 0;
}


