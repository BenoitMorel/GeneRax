#include "FamiliesFileParser.hpp"
#include <fstream>
#include <algorithm>
#include <IO/Logger.hpp>
#include <parallelization/ParallelContext.hpp>

enum FFPStep {
  header,
  reading_family
};


static bool update_family(const std::string &line, 
    FamilyInfo &currentFamily,
    Families &families)
{
  if (line[0] == '-') {
    if (currentFamily.name.size()) {
      families.push_back(currentFamily);
      currentFamily.reset();
    }
    currentFamily.name = line.substr(1, line.size() - 1);
    return true;
  }
  size_t delim_pos = line.find_first_of("=");
  std::string key = line.substr(0, delim_pos);
  std::string value = line.substr(delim_pos + 1, line.size() - delim_pos);
  if (key == "alignment") {
    currentFamily.alignmentFile = value;
  } else if (key == "starting_gene_tree") {
    currentFamily.startingGeneTree = value;
  } else if (key == "mapping") {
    currentFamily.mappingFile = value;
  } else if (key == "subst_model") {
    currentFamily.libpllModel = value;
  } else {
    Logger::error << "Unknown prefix " << key << std::endl;
    return false;
  }
  return true;
}

Families FamiliesFileParser::parseFamiliesFile(const std::string &familiesFile)
{
  Families families;
  std::ifstream reader(familiesFile);
  std::string line;
  FFPStep step = header;
  FamilyInfo currentFamily;
  int lineNumber = -1;
  while (getline(reader, line))  {
    lineNumber++;
    line.erase(remove_if(line.begin(), line.end(), ::isspace), line.end());
    line = line.substr(0, line.find("#"));
    if (!line.size() || line[0] == '#') {
      continue;
    }
    switch(step) {
    case header:
      if (line == "[FAMILIES]") {
        step = reading_family;
      }
      break;
    case reading_family:
      if (!update_family(line, currentFamily, families)) {
        Logger::error << "Error when parsing " << familiesFile << ":" << lineNumber << std::endl;
        ParallelContext::abort(1);
      }
      break;
    }
  }
  
  // terminate 
  if (step == reading_family) {
    families.push_back(currentFamily);
  }
  return families;
}

