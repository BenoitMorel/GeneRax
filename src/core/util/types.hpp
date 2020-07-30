#pragma once
#include <vector>
#include <unordered_map>

typedef std::vector<double> VectorDouble;
typedef std::vector<VectorDouble> MatrixDouble;
using VectorUint = std::vector<unsigned int>;
using MatrixUint = std::vector<VectorUint>;
typedef std::unordered_map<std::string, unsigned int> StringToUint;

