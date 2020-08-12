#include <IO/CustomNewickParser.hpp>
#include <string>
#include <iostream>
#include <chrono>
#include <vector>
#include <random>
#include <functional> 
#include <algorithm>  
#include <trees/PLLRootedTree.hpp>
#include <fstream>


void test_aux(const std::string &newickString, bool as_file)
{
  RTreeParsingError error;
  std::string input = newickString;
  if (as_file) {
    input = "temp.txt";
    std::ofstream os(input);
    os << newickString;
  }
  auto pllTree = as_file ? 
    pll_rtree_parse_newick(input.c_str()) 
    : pll_rtree_parse_newick_string(input.c_str());
  auto customTree = custom_rtree_parse_newick(input.c_str(), 
      as_file,
      &error);
  auto pllTreeString = pll_rtree_export_newick(pllTree->root, nullptr);
  auto customTreeString = pll_rtree_export_newick(customTree->root, nullptr);
  assert(std::string(pllTreeString) == std::string(customTreeString));
  free(pllTreeString);
  free(customTreeString);
  pll_rtree_destroy(pllTree, nullptr);
  pll_rtree_destroy(customTree, nullptr);
}

void test(const std::string &newickString)
{
  //test_aux(newickString, true);
  test_aux(newickString, false);
}

void benchmark_aux(const std::vector<std::string> &newickStrings,
    std::vector<pll_rtree_t *> &trees,
    bool pllVersion,
    bool as_file)
{
  for (auto &newickString: newickStrings) {
    pll_rtree_t *tree = nullptr;
    if (pllVersion) {
      if (as_file) {
        tree = pll_rtree_parse_newick(newickString.c_str());
      } else {
        tree = pll_rtree_parse_newick_string(newickString.c_str());
      }
    } else {
      RTreeParsingError error;
      tree = custom_rtree_parse_newick(newickString.c_str(), 
          as_file, 
          &error);
    }
    trees.push_back(tree);
  }
}

std::string random_string( size_t length )
{
  auto randchar = []() -> char
  {
    const char charset[] =
      "0123456789"
      "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
      "abcdefghijklmnopqrstuvwxyz";
    const size_t max_index = (sizeof(charset) - 1);
    return charset[ rand() % max_index ];
  };
  std::string str(length,0);
  std::generate_n( str.begin(), length, randchar );
  return str;
}

std::string getRandomNewick()
{
  std::unordered_set<std::string> labels;
  unsigned int labelsNumber =(rand() % 20000) + 4;
  unsigned int labelsSize = 5 + (rand() % 20);
  for (unsigned int i = 0; i < labelsNumber; ++i) {
    labels.insert(random_string(labelsSize));
  }
  PLLRootedTree tree(labels);
  return tree.getNewickString();
}
    

void benchmark()
{
  std::cerr << "[benchmark] generating random newick strings..." << std::endl;
  std::vector<std::string> newickStrings;
  unsigned int seed = 42;
  unsigned int trees = 200;
  unsigned int iterations = 100;
  srand(seed);
  for (unsigned int i = 0; i < trees; ++i) {
    newickStrings.push_back(getRandomNewick());
  }
  std::vector<pll_rtree_t *> pllTrees;
  std::vector<pll_rtree_t *> customTrees;
  auto startCustom = std::chrono::high_resolution_clock::now();
  std::cerr << "[benchmark] parsing with custom..." << std::endl;
  for (unsigned int i = 0; i < iterations; ++i) {
    benchmark_aux(newickStrings, customTrees, false, false);
  }
  auto elapsedCustom = std::chrono::high_resolution_clock::now() - startCustom;
  auto startPLL = std::chrono::high_resolution_clock::now();
  std::cerr << "[benchmark] parsing with pll..." << std::endl;
  for (unsigned int i = 0; i < iterations; ++i) {
    benchmark_aux(newickStrings, pllTrees, true, false);
  }
  auto elapsedPLL = std::chrono::high_resolution_clock::now() - startPLL;
  
  auto timePLL = std::chrono::duration_cast<std::chrono::seconds>(
      elapsedPLL).count();
  auto timeCustom = std::chrono::duration_cast<std::chrono::seconds>(
      elapsedCustom).count();
  std::cerr << "PLL parser: " << timePLL << "s" << std::endl;
  std::cerr << "Custom parser: " << timeCustom << "s" << std::endl;
  std::cerr << "Destroying all trees..." << std::endl;
  for (auto tree: pllTrees) {
    pll_rtree_destroy(tree, nullptr);
  }
  for (auto tree: customTrees) {
    pll_rtree_destroy(tree, nullptr);
  }
}


int main(int, char**)
{
   
  std::string easy = "((a,b)ab,(c,d)cd)root;";
  test(easy);
  std::string noInternalLabels = "((a,b),(c,d));";
  test(noInternalLabels);
  std::string numericLabels = "((a,b)13,(c,d)4cd)root;";
  test(numericLabels);
  std::string branchLengths = "((a:30.5,b:0.03)ab,(c:0,d:3)cd)root;";
  test(branchLengths);
  std::string branchLengthsScientific 
    = "((a:1e-10,b:0.03)ab,(c:0,d:3)cd)root;";
  test(branchLengthsScientific);
  benchmark();
  return 0;


}


