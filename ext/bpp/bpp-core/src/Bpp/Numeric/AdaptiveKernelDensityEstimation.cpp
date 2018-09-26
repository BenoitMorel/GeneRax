//
// File: AdaptiveKernelDensityEstimation.cpp
// Created by: Julien Dutheil
// Created on: Thu Nov 05 13:25:07 2009
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for numerical calculus. This file is part of the Bio++ project.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#include "AdaptiveKernelDensityEstimation.h"
#include "Matrix/MatrixTools.h"
#include "NumConstants.h"

using namespace bpp;
using namespace std;

void AdaptiveKernelDensityEstimation::init_()
{
  //Compute the covariance matrix of the sample:
  MatrixTools::covar(x_, covar_);

  //Compute the mean vector
  sampleMean_(x_, xMean_);
  
  //Compute the inverse of the square root of the covariance matrix:
  MatrixTools::pow<double>(covar_, -0.5, invSqrtCovar_);   

  //Compute the bandwidth:
  h_ = std::pow(4. / ((2 * static_cast<double>(r_) + 1.) * static_cast<double>(n_)), 1. / (static_cast<double>(r_) + 4.));
  //Compute as much as we can in advance to simplify the density calculation:
  c1_ = 1. / (std::sqrt(MatrixTools::det(covar_)) * static_cast<double>(n_) * std::pow(h_, static_cast<int>(r_)));
  
  //Now compute the local tuning of the bandwidth.
  //First estimate the pilot density:
  vector<double> xi(r_);
  LinearMatrix<double> diff_xi(r_, 1);
  LinearMatrix<double>  std_xi(r_, 1);
  for (unsigned int i = 0; i < n_; i++)
  {
    //Get the current xi point to evaluate:
    for(unsigned int k = 0; k < r_; k++)
      xi[k] = x_(k, i);
     
    //Sum loop, over all xi's:
    double sum = 0;
    for (unsigned int j = 0; j < n_; j++)
    {
      for (unsigned int k = 0; k < r_; k++)
        diff_xi(k, 0) = xi[k] - x_(k, j);
      MatrixTools::mult(invSqrtCovar_, diff_xi, std_xi);
      MatrixTools::scale(std_xi, 1. / h_);
      sum += kernel_(std_xi);
    }
    pilot_[i] = c1_ * sum;
  }

  //Compute the tuning parameters:
  double g = 0;
  for (unsigned int i = 0; i < n_; i++)
    g += std::log(pilot_[i]);
  g = std::exp(g / static_cast<double>(n_));
  for (unsigned int i = 0; i < n_; i++)
    lambda_[i] = std::pow(g / pilot_[i], gamma_);

  //Compute as much as we can in advance to simplify the density calculation:
  for (unsigned int i = 0; i < n_; i++)
    c2_[i] = std::pow(lambda_[i], - static_cast<double>(r_));
}

void AdaptiveKernelDensityEstimation::sampleMean_(const Matrix<double>& x, std::vector<double>& mean)
{
  size_t nc = x.getNumberOfColumns();
  size_t nr = x.getNumberOfRows();
  mean.resize(nr);
  for (size_t i = 0; i < nr; i++)
  {
    mean[i] = 0;
    for (size_t j = 0; j < nc; j++)
      mean[i] += x(i, j);
    mean[i] /= static_cast<double>(nc);
  }
}

double AdaptiveKernelDensityEstimation::kernel_(const Matrix<double>& x)
{
  //x is supposed to have only one column and r_ rows.
  //We compute the scalar product of the column with itself:
  double scalar = 0;
  for (size_t i = 0; i < r_; i++)
    scalar += std::pow(x(i, 0), 2.);

  return std::pow(2. * NumConstants::PI(), -static_cast<double>(r_) / 2.) * std::exp(-0.5 * scalar);
}

double AdaptiveKernelDensityEstimation::kDensity(const std::vector<double>& x)
{
  LinearMatrix<double> diff_xi(r_, 1);
  LinearMatrix<double> std_xi(r_, 1);
  //Sum loop, over all xi's:
  double sum = 0;
  for(unsigned int j = 0; j < n_; j++)
  {
    for(unsigned int k = 0; k < r_; k++)
      diff_xi(k, 0) = x[k] - x_(k, j);
    MatrixTools::mult(invSqrtCovar_, diff_xi, std_xi);
    MatrixTools::scale(std_xi, 1. / (h_ * lambda_[j]));
    sum += kernel_(std_xi) * c2_[j];
  }
  return c1_ * sum;  
}

