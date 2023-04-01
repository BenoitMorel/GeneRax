#include "HighwayCandidateParser.hpp"
#include <sstream>
#include <algorithm>
#include <IO/IO.hpp>


bool readTaxa(std::stringstream &iss, 
    std::vector<corax_rnode_t *> &nodes,
    const std::string &candidateFile,
    unsigned int lineNumber,
    const std::unordered_map<std::string, corax_rnode_t *> &labelToNode,
    bool from) 
{
  std::string str;
  char delimiter = from ? ',' : '\n';
  std::string fromToString = from ? "from" : "to";
  if (!std::getline(iss, str, delimiter)) {
    Logger::info << "Error: can't read the " << fromToString << " from species at line " << lineNumber << " of file " << candidateFile << std::endl;
    return false;
  }  
  IO::removeSpaces(str);
  if (str == "*") {
    for (auto it: labelToNode) {
      nodes.push_back(it.second);
    }
    return true;
  }
  auto it = labelToNode.find(str);
  if (it == labelToNode.end()) {
    Logger::info << "Error: " << fromToString << " label " << str << 
      " not found in the species tree. Check line " 
      << lineNumber << " of file " << candidateFile << std::endl;
    return false;
  }
  nodes.push_back(it->second);
  return true;
}

std::vector<Highway> HighwayCandidateParser::parse(
    const std::string &candidateFile,
    PLLRootedTree &speciesTree)
{
  std::vector<Highway> candidates;
  auto labelToNode = speciesTree.getLabelToNode(false);
  std::ifstream is(candidateFile);
  std::string line;
  size_t lineNumber = 0;
  while (std::getline(is, line)) {
    lineNumber++;
    if (IO::isBlanck(line)) {
      continue;
    }
    IO::removeSpaces(line);
    if (line[0] =='#') {
      continue;
    }
    bool ok = true;
    std::stringstream iss(line);
    std::vector<corax_rnode_t *> froms;
    std::vector<corax_rnode_t *> tos;
    ok &= readTaxa(iss, froms, candidateFile, lineNumber, labelToNode, true);
  
    ok &= readTaxa(iss, tos, candidateFile, lineNumber, labelToNode, false);
    if (!ok) {
      candidates.clear();
      return candidates;
    }
    for (auto from: froms) {
      for(auto to: tos) {
        candidates.push_back(Highway(from, to));
      }
    } 
  }

  return candidates;
}

