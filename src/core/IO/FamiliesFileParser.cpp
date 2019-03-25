#include "FamiliesFileParser.hpp"
#include <fstream>
#include <algorithm>
#include <IO/Logger.hpp>
#include <ParallelContext.hpp>

enum FFPStep {
  header,
  reading_family
};


bool update_family(const string line, 
    FamiliesFileParser::FamilyInfo &currentFamily,
    vector<FamiliesFileParser::FamilyInfo> &families)
{
  switch(line[0]) {
    case '-':
      if (currentFamily.name.size()) {
        families.push_back(currentFamily);
        currentFamily.reset();
      }
      currentFamily.name = line.substr(1, line.size() - 1);
      break;
    case 'A': 
      currentFamily.alignmentFile = line.substr(2, line.size() - 1);
      break;
    case 'G': 
      currentFamily.startingGeneTree = line.substr(2, line.size() - 1);
      break;
    case'M':
      currentFamily.mappingFile = line.substr(2, line.size() - 1);
      break;
    case'L':
      currentFamily.libpllModel = line.substr(2, line.size() - 1);
      Logger::info << "Libpll model: " << endl;
      break;
    default:
      return false;
  }
  return true;
}

vector<FamiliesFileParser::FamilyInfo> FamiliesFileParser::parseFamiliesFile(const string &familiesFile)
{
  vector<FamilyInfo> families;
  ifstream reader(familiesFile);
  string line;
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
        Logger::error << "Error when parsing " << familiesFile << ":" << lineNumber << endl;
        Logger::error << "The line should start with one of: - A G M" << endl;
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

