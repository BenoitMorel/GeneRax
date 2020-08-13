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
#include <unordered_map>

void test_aux(const std::string &newickString, bool as_file)
{
  RTreeParsingError error;
  std::string input = newickString;
  if (as_file) {
    input = "temp.txt";
    std::ofstream os(input);
    os << newickString;
    os.close();
  }
  auto pllTree = as_file ? 
    pll_rtree_parse_newick(input.c_str()) 
    : pll_rtree_parse_newick_string(input.c_str());
  auto customTree = custom_rtree_parse_newick(input.c_str(), 
      as_file,
      &error);
  assert(customTree);
  assert(pllTree);
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
  test_aux(newickString, true);
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
  unsigned int labelsNumber =(rand() % 2000) + 4;
  unsigned int labelsSize = 5 + (rand() % 20);
  for (unsigned int i = 0; i < labelsNumber; ++i) {
    labels.insert(random_string(labelsSize));
  }
  PLLRootedTree tree(labels);
  return tree.getNewickString();
}
    

void benchmark(bool asFile)
{
  if (asFile) {
    std::cerr << "Benchmarking from file" << std::endl;
  } else {
    std::cerr << "Benchmarking from string" << std::endl;
  }
  std::cerr << "[benchmark] generating random newick strings..." << std::endl;
  std::vector<std::string> newickStrings;
  unsigned int seed = 42;
  unsigned int trees = 100;
  unsigned int iterations = 100;
  srand(seed);
  for (unsigned int i = 0; i < trees; ++i) {
    auto newick = getRandomNewick();
    if (asFile) {
      std::string file("temp_");
      file += std::to_string(i);
      file += ".txt";
      std::ofstream os(file);
      os << newick;
      newickStrings.push_back(file);
    } else {
      newickStrings.push_back(newick);
    }
  }
  std::vector<pll_rtree_t *> pllTrees;
  std::vector<pll_rtree_t *> customTrees;
  auto startCustom = std::chrono::high_resolution_clock::now();
  std::cerr << "[benchmark] parsing with custom..." << std::endl;
  for (unsigned int i = 0; i < iterations; ++i) {
    benchmark_aux(newickStrings, customTrees, false, asFile);
  }
  auto elapsedCustom = std::chrono::high_resolution_clock::now() - startCustom;
  auto startPLL = std::chrono::high_resolution_clock::now();
  std::cerr << "[benchmark] parsing with pll..." << std::endl;
  for (unsigned int i = 0; i < iterations; ++i) {
    benchmark_aux(newickStrings, pllTrees, true, asFile);
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

void test_bad_trees_aux(const std::string &name,
    const std::string &tree, 
    ParsingErrorType expectedError)
{
  std::cout << "[Test bad tree] " << name << ": " << tree << std::endl;;
  RTreeParsingError error;
  bool as_file = false;
  auto customTree = custom_rtree_parse_newick(tree.c_str(), 
      as_file,
      &error);
  assert(error.type == expectedError);
  std::cout << "[Test bad tree]   OK!" << std::endl;
}

void test_good_trees()
{
  std::unordered_map<std::string, std::string> newicks;
  newicks.insert({"Easy", "((a,b)ab,(c,d)cd)root;"});
  newicks.insert({"No internal label", "((a,b),(c,d));"});
  newicks.insert({"Numeric labels", "((a,b)13,(c,d)4cd)root;"});
  newicks.insert({"Branch lengths", "((a:30.5,b:0.03)ab,(c:0,d:3)cd)root;"});
  newicks.insert({"Scientific notation", "((a:1e-10,b:0.03)ab,(c:0,d:3)cd)root;"});
  newicks.insert({"Spaces","( (a : 30.5 , b : 0.03 ) ab , (c :0,d : 3 ) cd )root;"});  
  newicks.insert({"Tabs","(\t(a\t:\t30.5\t,\tb\t:\t0.03\t)\tab\t,\t(c\t:0,d\t:\t3\t)\tcd\t)root\t;"});  
  
  for (auto &entry: newicks) {
    std::cout << "[Test good tree] " << entry.first << 
      ": " << entry.second << std::endl;
    test(entry.second);
    std::cout << "[Test good tree]   OK" << std::endl;
  }
}

void test_bad_trees()
{
  
  test_bad_trees_aux("Unrooted",
      "((a,b),(c,d),(e,f));", 
      PET_POLYTOMY);
  test_bad_trees_aux("Polytomy",
      "((a,b),(c,(d, e, f));", 
      PET_POLYTOMY);
  test_bad_trees_aux("No semicolon",
      "((a,b),(c,(d, e)))", 
      PET_NOSEMICOLON);
  test_bad_trees_aux("Invalid branch length",
      "((a,b),(c,(d, e:0.hi)));", 
      PET_INVALID_BRANCH_LENGTH);
  /*
  test_bad_trees_aux("Too many left parenthesis",
      "((a,b),(c,(d, e:0.1));", 
      PET_DOUBLE_BRANCH_LENGTH);
  test_bad_trees_aux("Too many right parenthesis",
      "(a,b),(c,(d, e:0.1)));", 
      PET_DOUBLE_BRANCH_LENGTH);
  test_bad_trees_aux("Double branch length",
      "((a,b),(c,(d, e:0.0:0.1)));", 
      PET_DOUBLE_BRANCH_LENGTH);
  */
}

int main(int, char**)
{
  test_good_trees();  
  test_bad_trees();
  //benchmark(false);
  //benchmark(true);
  return 0;


}


