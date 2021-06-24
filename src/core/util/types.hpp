#pragma once
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <array>

using VectorDouble = std::vector<double>;
using MatrixDouble = std::vector<VectorDouble>;
using DistanceMatrix = MatrixDouble;
using VectorUint = std::vector<unsigned int>;
using MatrixUint = std::vector<VectorUint>;
using StringToUint = std::unordered_map<std::string, unsigned int>;
using StringToInt = std::unordered_map<std::string, int>;

using SPID = unsigned int;
using TaxaSet = std::unordered_set<SPID>;
using MetaQuartet = std::array<TaxaSet, 4>;
using UInt3 = std::array<unsigned int, 3>;
using UInt4 = std::array<unsigned int, 4>;

using PerFamLL = std::vector<double>;
struct TransferFrequencies {
  MatrixUint count;
  std::vector<std::string> idToLabel;
};

