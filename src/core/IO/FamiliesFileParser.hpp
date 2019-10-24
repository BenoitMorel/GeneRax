#pragma once

#include <string>
#include <vector>


#include <families/Families.hpp>

class FamiliesFileParser {
public:
  FamiliesFileParser() = delete;
  static Families parseFamiliesFile(const std::string &familiesFile);
};
