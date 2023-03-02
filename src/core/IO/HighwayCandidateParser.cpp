#include "HighwayCandidateParser.hpp"
#include <sstream>
#include <algorithm>

static bool isBlanck(const std::string &s) 
{
  return s.empty() 
    || std::all_of(s.begin(), 
        s.end(), 
        [](char c){return std::isspace(c);});
}

template<typename T, typename P>
T remove_if(T beg, T end, P pred)
{
  T dest = beg;
  for (T itr = beg;itr != end; ++itr)
    if (!pred(*itr))
      *(dest++) = *itr;
  return dest;
}

void removeSpaces(std::string &str) {
  std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
  str.erase(end_pos, str.end());
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
    if (isBlanck(line)) {
      continue;
    }
    removeSpaces(line);
    if (line[0] =='#') {
      continue;
    }
    bool ok = true;
    std::string from;
    std::string to;
    std::stringstream iss(line);
    if (!std::getline(iss, from, ',')) {
      Logger::info << "Error: can't read the from species at line " << lineNumber << " of file " << candidateFile << std::endl;
      ok = false;
    }   
    if (ok && !std::getline(iss, to)) {
      Logger::info << "Error: can't read the dest species at line " << lineNumber << " of file " << candidateFile << std::endl;
      ok = false;
    }
    if (ok && labelToNode.find(from) == labelToNode.end()) {
      ok = false;
      Logger::info << "Error: from label " << from << 
        " not found in the species tree. Check your highway candidate file at line " 
        << lineNumber << " of file " << candidateFile << std::endl;
    }
    //removeSpaces(from);
    //removeSpaces(to);
    if (ok && labelToNode.find(to) == labelToNode.end()) {
      ok = false;
      Logger::info << "Error: to label " << to << 
        " not found in the species tree. Check your highway candidate file at line " 
        << lineNumber << " of file " << candidateFile << std::endl;
    }
    candidates.push_back(Highway(labelToNode[from], labelToNode[to]));
    if (!ok) {
      candidates.clear();
      return candidates;
    }
  }

  return candidates;
}

