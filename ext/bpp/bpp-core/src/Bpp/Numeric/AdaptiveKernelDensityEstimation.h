//
// File: AdaptiveKernelDensityEstimation.h
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

#ifndef _ADAPTIVEKERNELDENSITYESTIMATION_H_
#define _ADAPTIVEKERNELDENSITYESTIMATION_H_

#include "Matrix/Matrix.h"

namespace bpp
{

/**
 * @brief Density estimation using the adaptive kernel method.
 *
 * For now this implementation is quite restricted, more options may be implemented later...
 *
 * The source for this method can be found is the appendix of the following paper:
 * Ivan Kojadinovic, _Computational Statistics and Data Analaysis_ (2004), 46:269-294 
 *
 * @author Julien Dutheil
 */
class AdaptiveKernelDensityEstimation
{
  private:
    RowMatrix<double> x_; //The original sample
    size_t n_;
    size_t r_;
    RowMatrix<double> covar_; //The covariance matrix, used for the linear transformation
    RowMatrix<double> invSqrtCovar_; //The inverse of the square root of the covariance matrix, used for the linear transformation
    std::vector<double> xMean_;
    double gamma_; //Tune the effect of the pilot density.
    double c1_;
    std::vector<double> c2_;
    double h_; //The bandwidth.
    std::vector<double> lambda_; //The local tuning coefficient of the bandwidth.
    std::vector<double> pilot_; //The pilot density

  public:
    /**
     * @brief Build a new AdaptiveKernelDensityEstimation object.
     * @param x A mtrix contianing the sample point, one point per column.
     * The row of the matrix are the dimension of the sampled vectors, wich can be 1.
     * @param gamma Controls the influence of the pilot density. A value of 0
     * maximizes the impact of the pilot density, and hence corresponds to the standard
     * Kernel Density Estimation method. A value in ]0,1] allows a local adjustement of
     * the bandwith. The 0.5 value is commonly used.
     */
    AdaptiveKernelDensityEstimation(const Matrix<double>& x, double gamma = 0.5):
      x_(x), n_(x.getNumberOfColumns()), r_(x.getNumberOfRows()),
      covar_(), invSqrtCovar_(), xMean_(), gamma_(gamma),
      c1_(0), c2_(x.getNumberOfColumns()), h_(0),
      lambda_(x.getNumberOfColumns()), pilot_(x.getNumberOfColumns())
    {
      init_();
    }
    virtual ~AdaptiveKernelDensityEstimation() {}

  public:

    /**
     * @return The value of the estimated density for point x.
     * @param x The point where to estimate the density.
     */
    double kDensity(const std::vector<double>& x);

  private:
    void init_();

    void sampleMean_(const Matrix<double>& x, std::vector<double>& mean);
    
    /**
     * @brief The kernel function.
     *
     * For now a standard normal density is used, further options may be added later,
     * including the possibility to use your own function.
     *
     * @param x The point for which to compute the density, as a matrix with 1 column and r_ rows.
     * @return The value of the kernel function at the corresponding point.
     */
    double kernel_(const Matrix<double>& x);

};

} //End of namespace bpp.

#endif //_ADAPTIVEKERNELDENSITYESTIMATION_H_
