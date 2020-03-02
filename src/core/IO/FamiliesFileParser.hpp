#pragma once

#include <string>
#include <vector>


#include <IO/Families.hpp>

class FamiliesFileParser {
public:
  FamiliesFileParser() = delete;
  static Families parseFamiliesFile(const std::string &familiesFile);
};
