#include <getopt.h>
#include <iostream>
#include <ale/Constants.h>
#include <ale/tools/SpeciesGeneMapper.h>
#include <treeSearch/JointTree.h>
using namespace std;

inline void print_help() {

  cout << "Options:                                                                                              " << endl;
  cout << "   -h, --help" << endl;
  cout << "\tprint this help, then exit.\n" << endl;

  cout << "   -g, --genetree GENETREE_FILENAME" << endl;
  cout << "\tgene tree filename in Newick, NHX or PhyloXML format.\n" << endl;

  cout << "   -s, --speciestree SPECIESTREE_FILENAME" << endl;
  cout << "\tspecies tree filename in Newick, NHX or PhyloXML format.\n" << endl;

  cout << "   -a, --alignments ALIGNMENTS_FILENAME" << endl;
  cout << "\tfile containing the paths to the per-gene-family multiple alignments.\n" << endl;

  cout << "   -S, --smap SMAP_FILENAME" << endl;
  cout << "\tgene/species map filename." << endl;
  cout << "\tFirst column as the gene names and the second, species names.\n" << endl;

  cout << "   -d, --dupcost VALUE" << endl;
  cout << "\tgene duplication cost (default value = " << DEFAULT_DUPLICATION_COST <<").\n" << endl;

  cout << "   -l, --losscost VALUE" << endl;
  cout << "\tgene loss cost (default value = " << DEFAULT_LOSS_COST <<").\n" << endl;

  cout << "   -o, --output OUTPUT_FILENAME" << endl;
  cout << "\toutput filename (default: GENETREE_FILENAME_reconciliation.FORMAT).\n" << endl;
  
#if defined(_OPENMP) && !defined(__clang__)
    cout << "   -P, --parallelize" << endl;
    cout << "\tparallelize if it's possible.\n" << endl;
  #endif
}

struct Arguments {
  Arguments():
    dupCost(DEFAULT_DUPLICATION_COST),
    lossCost(DEFAULT_LOSS_COST),
    threads(1)
  {}

  string geneTree;
  string speciesTree;
  string alignments;
  string smap;
  MAPPING_METHOD mappingMethod;
  double dupCost;
  double lossCost;
  int threads;
  string output;
};

void parseArguments(int argc, char *argv[], Arguments &arguments) {

  static struct option long_options_lists[] = {
      {"help",              no_argument,       NULL,  'h'},
      {"genetree",          required_argument, NULL,  'g'},
      {"speciestree",       required_argument, NULL,  's'},
      {"alignments",        no_argument,       NULL,  'a'},
      {"smap",              required_argument, NULL,  'S'},
      {"dupcost",           required_argument, NULL,  'd'},
      {"losscost",          required_argument, NULL,  'l'},
      {"output",            required_argument, NULL,  'o'},
      #if defined(_OPENMP) && !defined(__clang__)
        {"parallelize",     no_argument,       NULL,  'P'},
      #endif
      {0,0,0,0}
  }; 
  int option;
  int option_index = 0;
  char options[50] = "hg:s:a:S:d:l:o:P:";
#if defined(_OPENMP) && !defined(__clang__)
    strcat(options, "P");
#endif
  
  if(argc == 1){
    print_help();
    exit(EXIT_SUCCESS);
  }

  while(((option = getopt_long(argc, argv, options, long_options_lists, &option_index))) != -1){
    switch(option){
      case 0:
      { // We are going to check if the long option given exists

        if (long_options_lists[option_index].flag != 0)
          break;
      }
      case 'h':
      {
        print_help();
        exit(EXIT_SUCCESS);
      }
      case 'g':
      {
        arguments.geneTree = string(optarg);
        break;
      }
      case 's':
      {
        arguments.speciesTree = string(optarg);
        break;
      }
      case 'a':
      {
        arguments.alignments = string(optarg);
        break;
      }
      case 'S':
      {
        arguments.smap = optarg;
        if (arguments.smap.find(".emf") != string::npos) {
          arguments.mappingMethod = ensembl;
        } else
          arguments.mappingMethod = smap;
        break;
      }
      case 'd':
      {
        arguments.dupCost = atof(optarg);
        break;
      }
      case 'l':
      {
        arguments.lossCost = atof(optarg);
        break;
      }
      case 'o':
      {
        arguments.output = optarg;
        break;
      }
      case 'P':
      {
        arguments.threads = atoi(optarg);
        break;
      }
      default:
      {
        cerr << "Syntax error" << endl;
        exit(EXIT_FAILURE);
      }
    }
  }

  if (arguments.geneTree == "") {
    cerr << "Please provide a gene tree" << endl;
    exit(EXIT_FAILURE);
  }
  if (arguments.speciesTree == "") {
    cerr << "Please provide a species tree" << endl;
    exit(EXIT_FAILURE);
  }
  if (arguments.alignments == "") {
    cerr << "Please provide an alignment file" << endl;
    exit(EXIT_FAILURE);
  }
  if (arguments.output == "") {
    cerr << "Please provide an output" << endl;
    exit(EXIT_FAILURE);
  }
}
  
void loadTrees(const string &filename, vector<BPPTree> &trees) {
  ifstream genefile(filename.c_str(), ios::in);
  string file_content((istreambuf_iterator<char>(genefile)), istreambuf_iterator<char>() ); // content of the file
  auto lines = Utils::splitString(file_content, "\n");
  for (auto &line: lines) {
    if (line[0] == '>') 
      continue;
    trees.push_back(BPPTree(IO::newickToPhyloTree(line, false)));
  }
}

int main(int argc, char * argv[]) {
  Arguments arguments;
  parseArguments(argc, argv, arguments);
  
  vector<BPPTree> geneTrees;
  vector<LibpllAlignmentInfo> alignmentsInfo;
  BPPTree speciesTree = IO::readTreeFile(arguments.speciesTree); 
  loadTrees(arguments.geneTree, geneTrees);
  auto maps = SpeciesGeneMapper::map(geneTrees.begin(), geneTrees.end(), 
       *speciesTree, trees, "",
       DEFAULT_SEPARATION_CHARACTER_BETWEEN_GENE_AND_SPECIES_NAME,
       DEFAULT_SPECIES_NAME_IS_PREFIX_IN_GENE_NAME,
       false);
  LibpllEvaluation::parseAlignmentInfo(arguments.alignments, alignmentsInfo, -1);
  if (alignmentsInfo.size() != geneTrees.size()) {
    cerr << "Error: the number of alignments (" << alignmentsInfo.size() << ") and gene trees (" << geneTrees.size() << ") differs" << endl;
  }
 
  /*
  ofstream ossizes (arguments.output + ".sizes");
  for (int i = 0; i < geneTrees.size(); ++i) {
    shared_ptr<LibpllEvaluation> eval = LibpllEvaluation::buildFromPhylo(geneTrees[i], alignmentsInfo[i]);
    auto partition = eval->getTreeInfo()->partitions[0];
    ossizes << partition->sites << " " << partition->tips << endl; 
  }
  */
  vector<double> likelihoods(geneTrees.size());
  #if defined(_OPENMP)
  # pragma omp parallel for schedule(dynamic) num_threads(arguments.threads)
  #endif
  for (int i = 0; i < geneTrees.size(); ++i) {
    JointTree jointTree(geneTrees[i], &(alignmentsInfo[i]), speciesTree, 
        maps[i], arguments.dupCost, arguments.lossCost);
    jointTree.optimizeParameters();
    double ll = jointTree.computeJointLoglk();
    #if defined(_OPENMP)
    # pragma omp critical
    #endif
    {
      likelihoods[i] = ll;
    }
  }
  ofstream os(arguments.output);
  double total_likelihood = 0.0;
  for (int i = 0; i < geneTrees.size(); ++i) {
      os << "tree" <<  i << " ll = " << likelihoods[i] << endl;
      total_likelihood += likelihoods[i];
  }
  os << "total trees ll = " << total_likelihood << endl;
  return 0;
}
