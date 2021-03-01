#pragma once

#include <vector>

using Likelihoods = std::vector<double>;

void performAUTest(std::vector<Likelihoods> likelihoods,
    size_t nboot,
    std::vector<double> &au_pvalues);


