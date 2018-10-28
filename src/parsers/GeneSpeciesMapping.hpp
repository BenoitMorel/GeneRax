#pragma once

#include <string>
#include <map>

using namespace std;

class GeneSpeciesMapping {
public:
  GeneSpeciesMapping(const string &mappingFile);
  const map<string, string> &getMap() const {return _map;}
  const string &getSpecies(const string &gene) {return _map[gene];}

private:
  map<string, string> _map;
};

