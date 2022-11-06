#include "LibpllParsers.hpp"
#include <fstream>
#include <streambuf>
#include <algorithm>
#include <parallelization/ParallelContext.hpp>
#include <cstring>
#include <sstream>
#include <stack>
#include <array>
#include <IO/Logger.hpp>
#include <IO/LibpllException.hpp>
#include <IO/RootedNewickParser.hpp>
#include <corax/corax.h>


static char * rtree_export_newick_recursive(const corax_rnode_t * root,
                                  char * (*cb_serialize)(const corax_rnode_t *))
{
  char * newick;
  assert(root != NULL);
  if (!(root->left) || !(root->right))
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
    }
    else
    {
      asprintf(&newick, "%s:%f", root->label, root->length);
    }
  }
  else
  {
    char * subtree1 = rtree_export_newick_recursive(root->left,cb_serialize);
    char * subtree2 = rtree_export_newick_recursive(root->right,cb_serialize);
    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      asprintf(&newick,
                              "(%s,%s)%s",
                              subtree1,
                              subtree2,
                              temp);
      free(temp);
    }
    else
    {
      asprintf(&newick,
                              "(%s,%s)%s:%f",
                              subtree1,
                              subtree2,
                              root->label ? root->label : "",
                              root->length);
    }
    free(subtree1);
    free(subtree2);
  }
  return newick;
}


char * corax_rtree_export_newick(const corax_rnode_t * root,
                                   char * (*cb_serialize)(const corax_rnode_t *))
{
  char * newick;
  if (!root) return NULL;

  if (!(root->left) || !(root->right))
  {
    if (cb_serialize)
    {
      newick = cb_serialize(root);
    }
    else
    {
      asprintf(&newick, "%s:%f", root->label, root->length);
    }
  }
  else
  {
    char * subtree1 = rtree_export_newick_recursive(root->left,cb_serialize);
    char * subtree2 = rtree_export_newick_recursive(root->right,cb_serialize);
    if (cb_serialize)
    {
      char * temp = cb_serialize(root);
      asprintf(&newick,
                              "(%s,%s)%s",
                              subtree1,
                              subtree2,
                              temp);
      free(temp);
    }
    else
    {
      asprintf(&newick,
                              "(%s,%s)%s:%f;",
                              subtree1,
                              subtree2,
                              root->label ? root->label : "",
                              root->length);
    }
    free(subtree1);
    free(subtree2);
  }

  return newick;
}

corax_utree_t *LibpllParsers::readNewickFromFile(const std::string &newickFilename)
{
  std::ifstream is(newickFilename);
  if (!is)
    throw LibpllException("Could not load open newick file ", newickFilename);

  std::string line;
  if (!std::getline(is, line)) {
    throw LibpllException("Error while reading tree (file is empty) from file: ", newickFilename); 
  }
  corax_utree_t *res = nullptr;
  try {
    res = readNewickFromStr(line);
  } catch (...) {
    throw LibpllException("Error while reading tree from file: ", newickFilename);
  }
  return res;
}


corax_utree_t *LibpllParsers::readNewickFromStr(const std::string &newickString)
{
  auto utree =  corax_utree_parse_newick_string_unroot(newickString.c_str());
  if (!utree) 
    throw LibpllException("Error while reading tree from std::string: ", newickString);
  return utree;
}


static corax_rtree_t *readRooted(const std::string &newick, bool isFile)
{
  ParsingError error;
  auto tree = custom_rtree_parse_newick(newick.c_str(), 
      isFile,
      &error);
  if (error.type == PET_NOERROR) {
    return tree;
  } else {
    std::string errorMessage;
    if (isFile) {
      errorMessage = "Error while reading rooted tree from file " + newick + ".\n";
    } else {
      errorMessage = "Error while reading rooted tree from string " + newick + ".\n";
    }
    errorMessage += "Error name: ";
    errorMessage += std::string(get_parsing_error_name(error.type));
    errorMessage += ".\n";
    errorMessage += "Error help message: ";
    errorMessage += std::string(get_parsing_error_diagnostic(error.type));
    errorMessage += ".\n";
    errorMessage += "The parsing error was detected at character ";
    errorMessage += std::to_string(error.offset) + ".";
    throw LibpllException(errorMessage); 
  }
}

corax_rtree_t *LibpllParsers::readRootedFromFile(const std::string &newickFile)
{
  return readRooted(newickFile, true);
}

corax_rtree_t *LibpllParsers::readRootedFromStr(const std::string &newickString)
{
  return readRooted(newickString, false);
}
  
void LibpllParsers::saveUtree(const corax_unode_t *utree, 
  const std::string &fileName, 
  bool append)
{
  std::ofstream os(fileName, (append ? std::ofstream::app : std::ofstream::out));
  char *newick = corax_utree_export_newick_rooted(utree, 0);
  os << newick;
  os.close();
  free(newick);
}
void LibpllParsers::saveRtree(const corax_rnode_t *rtree, 
    const std::string &fileName)
{
  std::ofstream os(fileName, std::ofstream::out);
  char *newick = corax_rtree_export_newick(rtree, 0);
  os << newick;
  os.close();
  free(newick);
}
  

void LibpllParsers::getUnodeNewickString(const corax_unode_t *rnode, std::string &newick)
{
  char *newickStr = corax_utree_export_newick(rnode, 0);
  newick = std::string(newickStr);
  free(newickStr);
}

void rtreeHierarchicalStringAux(const corax_rnode_t *node, std::vector<bool> &lefts, std::ostringstream &os)
{
  if (!node) {
    return;
  }
  for (unsigned int i = 0; i < lefts.size(); ++i) {
    auto left = lefts[i];
    if (i == lefts.size() - 1) {
      os << "---";
    } else {
      if (left) {
        os << "|  ";
      } else {
        os << "   ";
      }
    }
  }
  os << (node->label ? node->label : "null") << std::endl;
  lefts.push_back(true);
  rtreeHierarchicalStringAux(node->left, lefts, os);
  lefts[lefts.size() - 1] = false;
  rtreeHierarchicalStringAux(node->right, lefts, os);
  lefts.pop_back();
}

void LibpllParsers::getRtreeHierarchicalString(const corax_rtree_t *rtree, std::string &newick)
{
  std::ostringstream os;
  std::vector<bool> lefts;
  rtreeHierarchicalStringAux(rtree->root, lefts, os);
  newick = os.str();
}

std::vector<unsigned int> LibpllParsers::parallelGetTreeSizes(const Families &families) 
{
  unsigned int treesNumber = static_cast<unsigned int>(families.size());
  std::vector<unsigned int> localTreeSizes((treesNumber - 1 ) / ParallelContext::getSize() + 1, 0);
  for (auto i = ParallelContext::getBegin(treesNumber); i < ParallelContext::getEnd(treesNumber); i ++) {
    corax_utree_t *tree = LibpllParsers::readNewickFromFile(families[i].startingGeneTree);
    unsigned int taxa = tree->tip_count;
    localTreeSizes[i - ParallelContext::getBegin(treesNumber)] = taxa;
    corax_utree_destroy(tree, 0);
  }
  std::vector<unsigned int> treeSizes;
  ParallelContext::concatenateUIntVectors(localTreeSizes, treeSizes);
  treeSizes.erase(remove(treeSizes.begin(), treeSizes.end(), 0), treeSizes.end());
  assert(treeSizes.size() == families.size());
  return treeSizes;
}
void LibpllParsers::fillLeavesFromUtree(corax_utree_t *utree, std::unordered_set<std::string> &leaves)
{
  for (unsigned int i = 0; i < utree->tip_count + utree->inner_count; ++i) {
    auto node = utree->nodes[i];
    if (!node->next) {
      leaves.insert(std::string(node->label));
    }
  }
}

void LibpllParsers::fillLeavesFromRtree(corax_rtree_t *rtree, std::unordered_set<std::string> &leaves)
{
  for (unsigned int i = 0; i < rtree->tip_count + rtree->inner_count; ++i) {
    auto node = rtree->nodes[i];
    if (!node->left) {
      leaves.insert(std::string(node->label));
    }
  }
}

void LibpllParsers::parseMSA(const std::string &alignmentFilename, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights)
{
  if (!std::ifstream(alignmentFilename.c_str()).good()) {
    throw LibpllException("Alignment file " + alignmentFilename + "does not exist");
  }
  try {
    parseFasta(alignmentFilename.c_str(),
        stateMap, sequences, weights);
  } catch (...) {
    parsePhylip(alignmentFilename.c_str(),
        stateMap, sequences,
        weights);
  }
}

void LibpllParsers::parseFasta(const char *fastaFile, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights)
{
  auto reader = corax_fasta_open(fastaFile, corax_map_fasta);
  if (!reader) {
    corax_fasta_close(reader);
    throw LibpllException("Cannot parse fasta file ", fastaFile);
  }
  char * head;
  long head_len;
  char *seq;
  long seq_len;
  long seqno;
  int length;
  while (corax_fasta_getnext(reader, &head, &head_len, &seq, &seq_len, &seqno)) {
    sequences.push_back(PLLSequencePtr(new PLLSequence(head, seq, static_cast<unsigned int>(seq_len))));
    length = static_cast<int>(seq_len);
  }
  unsigned int count = static_cast<unsigned int>(sequences.size());
  char** buffer = static_cast<char**>(malloc(static_cast<size_t>(count) * sizeof(char *)));
  assert(buffer);
  for (unsigned int i = 0; i < count; ++i) {
    buffer[i] = sequences[i]->seq;
  }
  weights = corax_compress_site_patterns(buffer, stateMap, static_cast<int>(count), &length);
  if (!weights) {
    corax_fasta_close(reader);
    free(buffer);
    throw LibpllException("Error while parsing fasta: cannot compress sites from ", fastaFile);
  }
  for (unsigned int i = 0; i < count; ++i) {
    sequences[i]->len = static_cast<unsigned int>(length);
  }
  free(buffer);
  corax_fasta_close(reader);
}
  
void LibpllParsers::parsePhylip(const char *phylipFile, 
    const corax_state_t *stateMap,
    PLLSequencePtrs &sequences,
    unsigned int *&weights)
{
  assert(phylipFile);
  assert(stateMap);
  std::unique_ptr<corax_phylip_t, void (*)(corax_phylip_t*)> reader(corax_phylip_open(phylipFile, corax_map_phylip),
      corax_phylip_close);
  if (!reader) {
    throw LibpllException("Error while opening phylip file ", phylipFile);
  }
  corax_msa_t *msa = nullptr;
  // todobenoit check memory leaks when using the std::exception trick
  try {
    msa = corax_phylip_parse_interleaved(reader.get());
    if (!msa) {
      throw LibpllException("failed to parse ", phylipFile);
    }
  } catch (...) {
    std::unique_ptr<corax_phylip_t, void(*)(corax_phylip_t*)> 
      reader2(corax_phylip_open(phylipFile, corax_map_phylip), corax_phylip_close);
    msa = corax_phylip_parse_sequential(reader2.get());
    if (!msa) {
      throw LibpllException("failed to parse ", phylipFile);
    }
  }
  weights = corax_compress_site_patterns(msa->sequence, stateMap, msa->count, &msa->length);
  if (!weights) 
    throw LibpllException("Error while parsing fasta: cannot compress sites");
  for (auto i = 0; i < msa->count; ++i) {
    PLLSequencePtr seq(new PLLSequence(msa->label[i], msa->sequence[i], static_cast<unsigned int>(msa->length)));
    sequences.push_back(std::move(seq));
    // avoid freeing these buffers with corax_msa_destroy
    msa->label[i] = nullptr;
    msa->sequence[i] = nullptr;
  }
  corax_msa_destroy(msa);
}
  
std::unique_ptr<Model> LibpllParsers::getModel(const std::string &modelStrOrFilename)
{
  std::string modelStr = modelStrOrFilename;
  std::ifstream f(modelStr);
  if (f.good()) {
    getline(f, modelStr);
    modelStr = modelStr.substr(0, modelStr.find(","));
  }
  return std::make_unique<Model>(modelStr);
}

bool LibpllParsers::fillLabelsFromAlignment(const std::string &alignmentFilename, 
    const std::string& modelStrOrFilename,  
    std::unordered_set<std::string> &leaves)
{
  auto model = getModel(modelStrOrFilename);
  PLLSequencePtrs sequences;
  unsigned int *patternWeights = nullptr;
  bool res = true;
  try { 
    LibpllParsers::parseMSA(alignmentFilename, model->charmap(), sequences, patternWeights);
  } catch (...) {
    res = false;
  }
  free(patternWeights);
  for (auto &sequence: sequences) {
    if (leaves.find(sequence->label) != leaves.end()) {
      // duplicate taxa!
      return false;
    }
    leaves.insert(sequence->label);
  }
  return res;
}
  
unsigned int LibpllParsers::getMSALength(const std::string &alignmentFilename,
      const std::string &modelStrOrFilename)
{
  auto model = getModel(modelStrOrFilename);
  PLLSequencePtrs sequences;
  unsigned int *patternWeights = nullptr;
  try { 
    LibpllParsers::parseMSA(alignmentFilename, model->charmap(), sequences, patternWeights);
  } catch (...) {
    return 0;
  }
  if (sequences.size() == 0) {
    return 0;
  }
  unsigned int length = 0;
  for (unsigned int i = 0; i < sequences[0]->len; ++i) {
    length += patternWeights[i];
  }
  free(patternWeights);
  return length;

}

double LibpllParsers::getMSAEntropy(const std::string &alignmentFilename,
      const std::string &modelStrOrFilename)
{
  auto model = getModel(modelStrOrFilename);
  PLLSequencePtrs sequences;
  unsigned int *patternWeights = nullptr;
  try { 
    LibpllParsers::parseMSA(alignmentFilename, model->charmap(), sequences, patternWeights);
  } catch (...) {
    return 0;
  }
  if (sequences.size() == 0) {
    return 0;
  }
  unsigned int nonGaps = 0;
  for (auto &sequence: sequences) {
    for (unsigned int i = 0; i < sequence->len; ++i) {
      if (sequence->seq[i] != '-') {
        nonGaps += patternWeights[i];
      }
    }
  }
  free(patternWeights);
  return double(nonGaps) / double(sequences.size());

}
  
bool LibpllParsers::areLabelsValid(std::unordered_set<std::string> &leaves)
{
  std::array<bool, 256> forbiddenCharacter{}; // all set to false
  forbiddenCharacter[';'] = true;
  forbiddenCharacter[')'] = true;
  forbiddenCharacter['('] = true;
  forbiddenCharacter['['] = true;
  forbiddenCharacter[']'] = true;
  forbiddenCharacter[','] = true;
  forbiddenCharacter[':'] = true;
  forbiddenCharacter[';'] = true;
  for (auto &label: leaves) {
    for (auto c: label) {
      if (forbiddenCharacter[c]) {
        std::cerr << "invalid label: " << label << std::endl;
        return false;
      }
    }
  }
  return true;
}
  
void LibpllParsers::writeSuperMatrixFasta(const SuperMatrix &superMatrix,
      const std::string &outputFile)
{
  std::ofstream os(outputFile);
  for (auto &p: superMatrix) {
    auto &label = p.first;
    auto &sequence = p.second;
    os << ">" << label << std::endl;
    os << sequence << std::endl;
  }
}

