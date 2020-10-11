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

