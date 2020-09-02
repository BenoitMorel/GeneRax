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
  ParsingError error;
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
  assert(pllTree);
  auto customTree = custom_rtree_parse_newick(input.c_str(), 
      as_file,
      &error);
  if (error.type != PET_NOERROR) {
    std::cout << get_parsing_error_name(error.type) << " " 
      << error.offset << std::endl;
  }
  assert(error.type == PET_NOERROR);
  assert(customTree);
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
      ParsingError error;
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
  ParsingError error;
  bool as_file = false;
  auto customTree = custom_rtree_parse_newick(tree.c_str(), 
      as_file,
      &error);
  if(error.type != expectedError) {
    std::cerr << get_parsing_error_name(error.type) << std::endl;
  }
  if (customTree) {
    auto customTreeString = pll_rtree_export_newick(customTree->root, nullptr);
    std::cout << customTreeString << std::endl;
    free(customTreeString);
  }
  assert(!customTree && error.type == expectedError);
  std::cout << "[Test bad tree]   OK!" << std::endl;
}

void test_good_trees()
{
  std::unordered_map<std::string, std::string> newicks;
  newicks.insert({"Easy", "((a,b)ab,(c,d)cd)root;"});
  newicks.insert({"No internal label", "((a,b),(c,d));"});
  newicks.insert({"Numeric labels", "((a,b)13,(c,d)4cd)root;"});
  newicks.insert({"Branch lengths", "((a:30.5,b:0.03):48.0,(c:0,d:3)cd)root;"});
  newicks.insert({"Scientific notation", "((a:1e-10,b:0.03)ab,(c:0,d:3E-5)cd)root;"});
  newicks.insert({"Spaces","( (a : 30.5 , b : 0.03 ) ab , (c :0,d : 3 ) cd )root;"});  
  newicks.insert({"Tabs","(\t(a\t:\t30.5\t,\tb\t:\t0.03\t)\tab\t,\t(c\t:0,d\t:\t3\t)\tcd\t)root\t;"});  
  newicks.insert({"Newlines", "((a,b)ab\n,(c,d\n)cd)root;"});
  newicks.insert({"Windows newlines", "((a,b)ab\r\n,(c,d\r\n)cd)root;"});
  newicks.insert({"Weird characters", "((!a+7=5,b^o&)ab,($$£*c,d/\\?!_-|)cd)ro#~ot;"});
  
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
      PET_UNROOTED);
  test_bad_trees_aux("Polytomy",
      "((a,b),(c,(d, e, f));", 
      PET_POLYTOMY);
  test_bad_trees_aux("No semicolon",
      "((a,b),(c,(d, e)))", 
      PET_NOSEMICOLON);
  test_bad_trees_aux("Early semicolon",
      "((a,b);,(c,(d, e)))", 
      PET_INVALID_PARENTHESES);
  test_bad_trees_aux("Comments",
      "((a[comment],b),(c,(d, e):0.5));", 
      PET_INVALID_LABEL);
  test_bad_trees_aux("Missing comma",
      "((a,b)(c,(d, e):0.5));", 
      PET_POLYTOMY);
  test_bad_trees_aux("Token after the newick string",
      "((a,b),(c,(d, e)));wtf", 
      PET_TOKEN_AFTER_SEMICOLON);
  test_bad_trees_aux("Too many left parenthesis",
      "((a,b),(c,(d, e:0.1));", 
      PET_INVALID_PARENTHESES);
  test_bad_trees_aux("Too many right parenthesis",
      "(a,b),(c,(d, e:0.1)));", 
      PET_INVALID_PARENTHESES);
  test_bad_trees_aux("Space in label",
      "((a,b),(c,(d, e)hello world);", 
      PET_INVALID_LABEL);
  test_bad_trees_aux("Windows newline in label",
      "((a,b),(c,(d, e)hello\r\nworld);", 
      PET_INVALID_LABEL);
  test_bad_trees_aux("Unprintable character",
      "((a,b),(c,(d, e)hello\x01world);", 
      PET_INVALID_LABEL);
  test_bad_trees_aux("Unprintable character",
      "((a,b),(c,(d, e)hello\x18world);", 
      PET_INVALID_LABEL);
  test_bad_trees_aux("Newline in label",
      "((a,b),(c,(d, e)hello\nworld);", 
      PET_INVALID_LABEL);
  test_bad_trees_aux("Tab in label",
      "((a,b),(c,(d, e)hello\tworld);", 
      PET_INVALID_LABEL);
  test_bad_trees_aux("Space in label",
      "((a,b),(c,(d, e)hello world);", 
      PET_INVALID_LABEL);
  test_bad_trees_aux("0 children in non-terminal",
      "((a,b),(c,(d, ())));", 
      PET_EMPTY_NODE);
  test_bad_trees_aux("only 1 child",
      "((a,b),(c,(d, (e))));", 
      PET_ONLY_ONE_CHILD);
  test_bad_trees_aux("BL before label",
      "((a,b),(c,(d, e):0.5 label));", 
      PET_INVALID_BRANCH_LENGTH);
  test_bad_trees_aux("Invalid BL",
      "((a,b),(c,(d, e):0.a5));", 
      PET_INVALID_BRANCH_LENGTH);
  test_bad_trees_aux("Double branch length",
      "((a,b),(c,(d, e:0.0:0.1)));", 
      PET_DOUBLE_BRANCH_LENGTH);
}

void produce_separator_table()
{
  std::cout << "{";
  std::string separators("()[]\"';:,\n\r \t");
  for (unsigned int i = 0; i < 256; ++i) {
    char c = char(i);
    bool isSeparator = false;
    isSeparator |= (i < 32);
    isSeparator |= separators.find(c) != std::string::npos;
    std::cout << isSeparator;
    if (i != 255) {
      std::cout << ",";
    }
  }
  std::cout << "};" << std::endl;
}

int main(int, char**)
{
  //produce_separator_table();
  //benchmark(false);
  //benchmark(true);
  test_good_trees();  
  test_bad_trees();
  return 0;


}


