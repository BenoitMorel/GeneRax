#ifndef _JOINTSEARCH_ARGUMENTS_HPP
#define _JOINTSEARCH_ARGUMENTS_HPP

#include <string>

using namespace std;

class Arguments {
public:
  static void init(int argc, char * argv[]);
  static void checkInputs();
  static void printHelp();
  static void printCommand();
  static void printSummary();
public:
  static int argc;
  static char ** argv;
  static string geneTree;
  static string alignment;
  static string speciesTree;
  static string strategy;
  static int threads;
  static string output;
  static bool check;
  static bool verbose;
};

#endif
