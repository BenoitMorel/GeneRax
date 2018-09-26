//
// File: EigenValue.h
// Created by: Julien Dutheil
// Created on: Tue Apr 7 16:24 2004
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



#ifndef _EIGENVALUE_H_
#define _EIGENVALUE_H_
#define TOST(i) static_cast<size_t>(i)

#include <algorithm>
// for min(), max() below

#include <cmath>
// for abs() below
#include <climits>

#include "Matrix.h"
#include "../NumTools.h"

namespace bpp
{

/** 
 * @brief Computes eigenvalues and eigenvectors of a real (non-complex) matrix. 
 * 
 * [This class and its documentation is adpated from the C++ port of the JAMA library.]
 * 
 * If A is symmetric, then A = V*D*V' where the eigenvalue matrix D is
 * diagonal and the eigenvector matrix V is orthogonal. That is,
 * the diagonal values of D are the eigenvalues, and
 * V*V' = I, where I is the identity matrix.  The columns of V 
 * represent the eigenvectors in the sense that A*V = V*D.
 *
 * If A is not symmetric, then the eigenvalue matrix D is block diagonal
 * with the real eigenvalues in 1-by-1 blocks and any complex eigenvalues,
 * a + i*b, in 2-by-2 blocks, [a, b; -b, a].  That is, if the complex
 * eigenvalues look like
 * <pre>
 * 
 *         u + iv     .        .          .      .    .
 *           .      u - iv     .          .      .    .
 *           .        .      a + ib       .      .    .
 *           .        .        .        a - ib   .    .
 *           .        .        .          .      x    .
 *           .        .        .          .      .    y
 * </pre>
 * then D looks like
 * <pre>
 * 
 *           u        v        .          .      .    .
 *          -v        u        .          .      .    . 
 *           .        .        a          b      .    .
 *           .        .       -b          a      .    .
 *           .        .        .          .      x    .
 *           .        .        .          .      .    y
 * </pre>
 * This keeps V a real matrix in both symmetric and non-symmetric
 *  cases, and A*V = V*D.
 *
 *
 * The matrix V may be badly
 * conditioned, or even singular, so the validity of the equation
 * A = V*D*inverse(V) depends upon the condition number of V.
 *
 * (Adapted from JAMA, a Java Matrix Library, developed by jointly 
 *  by the Mathworks and NIST; see  http://math.nist.gov/javanumerics/jama).
 */
template <class Real>
class EigenValue
{

  private:
   /** 
    * @brief Row and column dimension (square matrix).
    */
   size_t n_;

   /**
    * @brief Tell if the matrix is symmetric.
    */
   bool issymmetric_;

   /**
    * @name Arrays for internal storage of eigenvalues.
    *
    * @{
    */
   std::vector<Real> d_;         /* real part */
   std::vector<Real> e_;         /* img part */
   /** @} */
   
   /**
    * @brief Array for internal storage of eigenvectors.
    */
   RowMatrix<Real> V_;

   /**
    * @brief Matrix for internal storage of nonsymmetric Hessenberg form.
    *
    * Internal storage of nonsymmetric Hessenberg form.
    */
   RowMatrix<Real> H_;
   

   /**
    * @brief Matrix for internal storage of eigen values in a matrix form.
    *
    * Internal storage of eigen values in a matrix form.
    */
   mutable RowMatrix<Real> D_;

   /**
    * @brief Working storage for nonsymmetric algorithm.
    *
    * Working storage for nonsymmetric algorithm.
    */
   std::vector<Real> ort_;

   /**
    * @brief Symmetric Householder reduction to tridiagonal form.
    *
    * This is derived from the Algol procedures tred2 by
    * Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    * Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    * Fortran subroutine in EISPACK.
    */
   void tred2()
   {
     for (size_t j = 0; j < n_; j++)
     {
       d_[j] = V_(n_-1,j);
     }

     // Householder reduction to tridiagonal form.
   
     for (size_t i = n_-1; i > 0; i--)
     {
       // Scale to avoid under/overflow.
   
       Real scale = 0.0;
       Real h = 0.0;
       for (size_t k = 0; k < i; ++k)
       {
         scale = scale + NumTools::abs<Real>(d_[k]);
       }
       if (scale == 0.0)
       {
         e_[i] = d_[i-1];
         for (size_t j = 0; j < i; ++j)
         {
           d_[j] = V_(i - 1, j);
           V_(i, j) = 0.0;
           V_(j, i) = 0.0;
         }
       }
       else
       {
         // Generate Householder vector.
   
         for (size_t k = 0; k < i; ++k)
         {
           d_[k] /= scale;
           h += d_[k] * d_[k];
         }
         Real f = d_[i - 1];
         Real g = sqrt(h);
         if (f > 0)
         {
           g = -g;
         }
         e_[i] = scale * g;
         h = h - f * g;
         d_[i-1] = f - g;
         for (size_t j = 0; j < i; ++j)
         {
           e_[j] = 0.0;
         }
   
         // Apply similarity transformation to remaining columns.
   
         for (size_t j = 0; j < i; ++j)
         {
           f = d_[j];
           V_(j,i) = f;
           g = e_[j] + V_(j,j) * f;
           for (size_t k = j + 1; k <= i - 1; k++)
           {
             g += V_(k,j) * d_[k];
             e_[k] += V_(k,j) * f;
           }
           e_[j] = g;
         }
         f = 0.0;
         for (size_t j = 0; j < i; ++j)
         {
           e_[j] /= h;
           f += e_[j] * d_[j];
         }
         Real hh = f / (h + h);
         for (size_t j = 0; j < i; ++j)
         {
           e_[j] -= hh * d_[j];
         }
         for (size_t j = 0; j < i; ++j)
         {
           f = d_[j];
           g = e_[j];
           for (size_t k = j; k <= i-1; ++k)
           {
             V_(k,j) -= (f * e_[k] + g * d_[k]);
           }
           d_[j] = V_(i-1,j);
           V_(i,j) = 0.0;
         }
       }
       d_[i] = h;
     }
   
     // Accumulate transformations.
   
     for (size_t i = 0; i < n_-1; i++)
     {
       V_(n_-1,i) = V_(i,i);
       V_(i,i) = 1.0;
       Real h = d_[i+1];
       if (h != 0.0)
       {
         for (size_t k = 0; k <= i; k++)
         {
           d_[k] = V_(k,i+1) / h;
         }
         for (size_t j = 0; j <= i; j++)
         {
           Real g = 0.0;
           for (size_t k = 0; k <= i; k++)
           {
             g += V_(k,i+1) * V_(k,j);
           }
           for (size_t k = 0; k <= i; k++)
           {
             V_(k,j) -= g * d_[k];
           }
         }
       }
       for (size_t k = 0; k <= i; k++)
       {
         V_(k,i+1) = 0.0;
       }
     }
     for (size_t j = 0; j < n_; j++)
     {
       d_[j] = V_(n_-1,j);
       V_(n_-1,j) = 0.0;
     }
     V_(n_-1,n_-1) = 1.0;
     e_[0] = 0.0;
   } 

   /**
    * @brief Symmetric tridiagonal QL algorithm.
    *
    * This is derived from the Algol procedures tql2, by
    * Bowdler, Martin, Reinsch, and Wilkinson, Handbook for
    * Auto. Comp., Vol.ii-Linear Algebra, and the corresponding
    * Fortran subroutine in EISPACK.
    */
   void tql2 ()
   {
     for (size_t i = 1; i < n_; i++)
     {
       e_[i - 1] = e_[i];
     }
     e_[n_ - 1] = 0.0;
   
     Real f = 0.0;
     Real tst1 = 0.0;
     Real eps = pow(2.0,-52.0);
     for (size_t l = 0; l < n_; ++l)
     {
       // Find small subdiagonal element
   
       tst1 = std::max(tst1, NumTools::abs<Real>(d_[l]) + NumTools::abs<Real>(e_[l]));
       size_t m = l;

       // Original while-loop from Java code
       while (m < n_)
       {
         if (NumTools::abs<Real>(e_[m]) <= eps*tst1)
         {
           break;
         }
         m++;
       }

       // If m == l, d_[l] is an eigenvalue,
       // otherwise, iterate.
   
       if (m > l)
       {
         int iter = 0;
         do
         {
           iter = iter + 1;  // (Could check iteration count here.)
   
           // Compute implicit shift
   
           Real g = d_[l];
           Real p = (d_[l + 1] - g) / (2.0 * e_[l]);
           Real r = hypot(p,1.0);
           if (p < 0)
           {
             r = -r;
           }
           d_[l] = e_[l] / (p + r);
           d_[l + 1] = e_[l] * (p + r);
           Real dl1 = d_[l + 1];
           Real h = g - d_[l];
           for (size_t i = l + 2; i < n_; ++i)
           {
             d_[i] -= h;
           }
           f = f + h;
   
           // Implicit QL transformation.
   
           p = d_[m];
           Real c = 1.0;
           Real c2 = c;
           Real c3 = c;
           Real el1 = e_[l + 1];
           Real s = 0.0;
           Real s2 = 0.0;
           //for (size_t i = m - 1; i >= l; --i)
           for (size_t ii = m; ii > l; --ii)
           {
             size_t i = ii - 1; //to avoid infinite loop!
             c3 = c2;
             c2 = c;
             s2 = s;
             g = c * e_[i];
             h = c * p;
             r = hypot(p, e_[i]);
             e_[i + 1] = s * r;
             s = e_[i] / r;
             c = p / r;
             p = c * d_[i] - s * g;
             d_[i + 1] = h + s * (c * g + s * d_[i]);
   
             // Accumulate transformation.
   
             for (size_t k = 0; k < n_; k++)
             {
               h = V_(k, i + 1);
               V_(k, i + 1) = s * V_(k, i) + c * h;
               V_(k, i) = c * V_(k, i) - s * h;
             }
           }
           p = -s * s2 * c3 * el1 * e_[l] / dl1;
           e_[l] = s * p;
           d_[l] = c * p;
   
           // Check for convergence.
   
         } while (NumTools::abs<Real>(e_[l]) > eps * tst1);
       }
       d_[l] = d_[l] + f;
       e_[l] = 0.0;
     }
     
     // Sort eigenvalues and corresponding vectors.
   
     for (size_t i = 0; n_ > 0 && i < n_-1; i++)
     {
       size_t k = i;
       Real p = d_[i];
       for (size_t j = i+1; j < n_; j++)
       {
         if (d_[j] < p)
         {
           k = j;
           p = d_[j];
         }
       }
       if (k != i)
       {
         d_[k] = d_[i];
         d_[i] = p;
         for (size_t j = 0; j < n_; j++)
         {
           p = V_(j, i);
           V_(j,i) = V_(j, k);
           V_(j,k) = p;
         }
       }
     }
   }

   /**
    * @brief Nonsymmetric reduction to Hessenberg form.
    *
    * This is derived from the Algol procedures orthes and ortran,
    * by Martin and Wilkinson, Handbook for Auto. Comp.,
    * Vol.ii-Linear Algebra, and the corresponding
    * Fortran subroutines in EISPACK.
    */
   void orthes()
   {
     if (n_ == 0) return;
     size_t low = 0;
     size_t high = n_-1;
   
     for (size_t m = low + 1; m <= high - 1; ++m)
     {
       // Scale column.
   
       Real scale = 0.0;
       for (size_t i = m; i <= high; ++i)
       {
         scale = scale + NumTools::abs<Real>(H_(i, m - 1));
       }
       if (scale != 0.0)
       {
         // Compute Householder transformation.
   
         Real h = 0.0;
         for (size_t i = high; i >= m; --i)
         {
           ort_[i] = H_(i, m - 1)/scale;
           h += ort_[i] * ort_[i];
         }
         Real g = sqrt(h);
         if (ort_[m] > 0)
         {
           g = -g;
         }
         h = h - ort_[m] * g;
         ort_[m] = ort_[m] - g;
   
         // Apply Householder similarity transformation
         // H = (I-u*u'/h)*H*(I-u*u')/h)
   
         for (size_t j = m; j < n_; ++j)
         {
           Real f = 0.0;
           for (size_t i = high; i >= m; --i)
           {
             f += ort_[i] * H_(i, j);
           }
           f = f/h;
           for (size_t i = m; i <= high; ++i)
           {
             H_(i,j) -= f*ort_[i];
           }
         }
   
         for (size_t i = 0; i <= high; ++i)
         {
           Real f = 0.0;
           for (size_t j = high; j >= m; --j)
           {
             f += ort_[j] * H_(i, j);
           }
           f = f/h;
           for (size_t j = m; j <= high; ++j)
           {
             H_(i,j) -= f * ort_[j];
           }
         }
         ort_[m] = scale * ort_[m];
         H_(m, m - 1) = scale * g;
       }
     }
   
     // Accumulate transformations (Algol's ortran).

     for (size_t i = 0; i < n_; i++)
     {
       for (size_t j = 0; j < n_; j++)
       {
         V_(i,j) = (i == j ? 1.0 : 0.0);
       }
     }

     for (size_t m = high - 1; m >= low + 1; --m)
     {
       if (H_(m, m - 1) != 0.0)
       {
         for (size_t i = m + 1; i <= high; ++i)
         {
           ort_[i] = H_(i, m - 1);
         }
         for (size_t j = m; j <= high; ++j)
         {
           Real g = 0.0;
           for (size_t i = m; i <= high; i++)
           {
             g += ort_[i] * V_(i, j);
           }
           // Double division avoids possible underflow
           g = (g / ort_[m]) / H_(m, m - 1);
           for (size_t i = m; i <= high; ++i)
           {
             V_(i, j) += g * ort_[i];
           }
         }
       }
     }
   }


   // Complex scalar division.

   Real cdivr, cdivi;
   void cdiv(Real xr, Real xi, Real yr, Real yi)
   {
     Real r,d;
     if (NumTools::abs<Real>(yr) > NumTools::abs<Real>(yi))
     {
       r = yi/yr;
       d = yr + r*yi;
       cdivr = (xr + r*xi)/d;
       cdivi = (xi - r*xr)/d;
     }
     else
     {
       r = yr/yi;
       d = yi + r*yr;
       cdivr = (r*xr + xi)/d;
       cdivi = (r*xi - xr)/d;
     }
   }


   // Nonsymmetric reduction from Hessenberg to real Schur form.

  void hqr2 ()
  {
    //  This is derived from the Algol procedure hqr2,
    //  by Martin and Wilkinson, Handbook for Auto. Comp.,
    //  Vol.ii-Linear Algebra, and the corresponding
    //  Fortran subroutine in EISPACK.
  
    // Initialize
   
    int nn = static_cast<int>(this->n_);
    int n = nn-1;
    int low = 0;
    int high = nn-1;
    Real eps = pow(2.0,-52.0);
    Real exshift = 0.0;
    Real p=0,q=0,r=0,s=0,z=0,t,w,x,y;
   
    // Store roots isolated by balanc and compute matrix norm
   
    Real norm = 0.0;
    for (int i = 0; i < nn; i++)
    {
      if ((i < low) || (i > high))
      {
        d_[TOST(i)] = H_(TOST(i),TOST(i));
        e_[TOST(i)] = 0.0;
      }
      for (int j = std::max(i-1,0); j < nn; j++)
      {
        norm = norm + NumTools::abs<Real>(H_(TOST(i),TOST(j)));
      }
    }
   
    // Outer loop over eigenvalue index
  
    int iter = 0;
    while (n >= low)
    {
      // Look for single small sub-diagonal element
  
      int l = n;
      while (l > low)
      {
        s = NumTools::abs<Real>(H_(TOST(l-1),TOST(l-1))) + NumTools::abs<Real>(H_(TOST(l),TOST(l)));
        if (s == 0.0)
        {
          s = norm;
        }
        if (NumTools::abs<Real>(H_(TOST(l),TOST(l-1))) < eps * s)
        {
          break;
        }
        l--;
      }
      
      // Check for convergence
      // One root found
  
      if (l == n)
      {
        H_(TOST(n),TOST(n)) = H_(TOST(n),TOST(n)) + exshift;
        d_[TOST(n)] = H_(TOST(n),TOST(n));
        e_[TOST(n)] = 0.0;
        n--;
        iter = 0;
  
        // Two roots found
  
      }
      else if (l == n-1)
      {
        w = H_(TOST(n),TOST(n-1)) * H_(TOST(n-1),TOST(n));
        p = (H_(TOST(n-1),TOST(n-1)) - H_(TOST(n),TOST(n))) / 2.0;
        q = p * p + w;
        z = sqrt(NumTools::abs<Real>(q));
        H_(TOST(n),TOST(n)) = H_(TOST(n),TOST(n)) + exshift;
        H_(TOST(n-1),TOST(n-1)) = H_(TOST(n-1),TOST(n-1)) + exshift;
        x = H_(TOST(n),TOST(n));
  
        // Real pair
  
        if (q >= 0)
        {
          if (p >= 0)
          {
            z = p + z;
          }
          else
          {
            z = p - z;
          }
          d_[TOST(n - 1)] = x + z;
          d_[TOST(n)] = d_[TOST(n - 1)];
          if (z != 0.0)
          {
            d_[TOST(n)] = x - w / z;
          }
          e_[TOST(n - 1)] = 0.0;
          e_[TOST(n)] = 0.0;
          x = H_(TOST(n),TOST(n-1));
          s = NumTools::abs<Real>(x) + NumTools::abs<Real>(z);
          p = x / s;
          q = z / s;
          r = sqrt(p * p+q * q);
          p = p / r;
          q = q / r;
  
          // Row modification
  
          for (int j = n-1; j < nn; j++)
          {
            z = H_(TOST(n-1),TOST(j));
            H_(TOST(n-1),TOST(j)) = q * z + p * H_(TOST(n),TOST(j));
            H_(TOST(n),TOST(j)) = q * H_(TOST(n),TOST(j)) - p * z;
          }
  
          // Column modification
  
          for (int i = 0; i <= n; i++)
          {
            z = H_(TOST(i),TOST(n-1));
            H_(TOST(i),TOST(n-1)) = q * z + p * H_(TOST(i),TOST(n));
            H_(TOST(i),TOST(n)) = q * H_(TOST(i),TOST(n)) - p * z;
          }
  
          // Accumulate transformations
  
          for (int i = low; i <= high; i++)
          {
            z = V_(TOST(i),TOST(n-1));
            V_(TOST(i),TOST(n-1)) = q * z + p * V_(TOST(i),TOST(n));
            V_(TOST(i),TOST(n)) = q * V_(TOST(i),TOST(n)) - p * z;
          }
   
          // Complex pair
  
        }
        else
        {
          d_[TOST(n - 1)] = x + p;
          d_[TOST(n)] = x + p;
          e_[TOST(n-1)] = z;
          e_[TOST(n)] = -z;
        }
        n = n - 2;
        iter = 0;
  
        // No convergence yet
  
      }
      else
      {
        // Form shift
  
        x = H_(TOST(n),TOST(n));
        y = 0.0;
        w = 0.0;
        if (l < n)
        {
          y = H_(TOST(n-1),TOST(n-1));
          w = H_(TOST(n),TOST(n-1)) * H_(TOST(n-1),TOST(n));
        }
  
        // Wilkinson's original ad hoc shift
  
        if (iter == 10)
        {
          exshift += x;
          for (int i = low; i <= n; i++)
          {
            H_(TOST(i),TOST(i)) -= x;
          }
          s = NumTools::abs<Real>(H_(TOST(n),TOST(n-1))) + NumTools::abs<Real>(H_(TOST(n-1),TOST(n-2)));
          x = y = 0.75 * s;
          w = -0.4375 * s * s;
        }

        // MATLAB's new ad hoc shift
        if (iter == 30)
        {
          s = (y - x) / 2.0;
          s = s * s + w;
          if (s > 0)
          {
            s = sqrt(s);
            if (y < x)
            {
              s = -s;
            }
            s = x - w / ((y - x) / 2.0 + s);
            for (int i = low; i <= n; i++)
            {
              H_(TOST(i),TOST(i)) -= s;
            }
            exshift += s;
            x = y = w = 0.964;
          }
        }
  
        iter = iter + 1;   // (Could check iteration count here.)
  
        // Look for two consecutive small sub-diagonal elements
  
        int m = n-2;
        while (m >= l)
        {
          z = H_(TOST(m),TOST(m));
          r = x - z;
          s = y - z;
          p = (r * s - w) / H_(TOST(m+1),TOST(m)) + H_(TOST(m),TOST(m+1));
          q = H_(TOST(m+1),TOST(m+1)) - z - r - s;
          r = H_(TOST(m+2),TOST(m+1));
          s = NumTools::abs<Real>(p) + NumTools::abs<Real>(q) + NumTools::abs<Real>(r);
          p = p / s;
          q = q / s;
          r = r / s;
          if (m == l)
          {
            break;
          }
          if (NumTools::abs<Real>(H_(TOST(m),TOST(m-1))) * (NumTools::abs<Real>(q) + NumTools::abs<Real>(r)) <
               eps * (NumTools::abs<Real>(p) * (NumTools::abs<Real>(H_(TOST(m-1),TOST(m-1))) + NumTools::abs<Real>(z) +
               NumTools::abs<Real>(H_(TOST(m+1),TOST(m+1))))))
          {
            break;
          }
          m--;
        }
  
        for (int i = m+2; i <= n; i++)
        {
          H_(TOST(i),TOST(i-2)) = 0.0;
          if (i > m+2)
          {
            H_(TOST(i),TOST(i-3)) = 0.0;
          }
        }
 
        // Double QR step involving rows l:n and columns m:n
   
        for (int k = m; k <= n-1; k++)
        {
          int notlast = (k != n-1);
          if (k != m)
          {
            p = H_(TOST(k),TOST(k-1));
            q = H_(TOST(k+1),TOST(k-1));
            r = (notlast ? H_(TOST(k+2),TOST(k-1)) : 0.0);
            x = NumTools::abs<Real>(p) + NumTools::abs<Real>(q) + NumTools::abs<Real>(r);
            if (x != 0.0)
            {
              p = p / x;
              q = q / x;
              r = r / x;
            }
          }
          if (x == 0.0)
          {
            break;
          }
          s = sqrt(p * p + q * q + r * r);
          if (p < 0)
          {
            s = -s;
          }
          if (s != 0)
          {
            if (k != m)
            {
              H_(TOST(k),TOST(k-1)) = -s * x;
            }
            else if (l != m)
            {
              H_(TOST(k),TOST(k-1)) = -H_(TOST(k),TOST(k-1));
            }
            p = p + s;
            x = p / s;
            y = q / s;
            z = r / s;
            q = q / p;
            r = r / p;
   
            // Row modification
   
            for (int j = k; j < nn; j++)
            {
              p = H_(TOST(k),TOST(j)) + q * H_(TOST(k+1),TOST(j));
              if (notlast)
              {
                p = p + r * H_(TOST(k+2),TOST(j));
                H_(TOST(k+2),TOST(j)) = H_(TOST(k+2),TOST(j)) - p * z;
              }
              H_(TOST(k),TOST(j)) = H_(TOST(k),TOST(j)) - p * x;
              H_(TOST(k+1),TOST(j)) = H_(TOST(k+1),TOST(j)) - p * y;
            }
   
            // Column modification
   
            for (int i = 0; i <= std::min(n,k+3); i++)
            {
              p = x * H_(TOST(i),TOST(k)) + y * H_(TOST(i),TOST(k+1));
              if (notlast)
              {
                p = p + z * H_(TOST(i),TOST(k+2));
                H_(TOST(i),TOST(k+2)) = H_(TOST(i),TOST(k+2)) - p * r;
              }
              H_(TOST(i),TOST(k)) = H_(TOST(i),TOST(k)) - p;
              H_(TOST(i),TOST(k+1)) = H_(TOST(i),TOST(k+1)) - p * q;
            }
   
            // Accumulate transformations
   
            for (int i = low; i <= high; i++)
            {
              p = x * V_(TOST(i),TOST(k)) + y * V_(TOST(i),TOST(k+1));
              if (notlast)
              {
                p = p + z * V_(TOST(i),TOST(k+2));
                V_(TOST(i),TOST(k+2)) = V_(TOST(i),TOST(k+2)) - p * r;
              }
              V_(TOST(i),TOST(k)) = V_(TOST(i),TOST(k)) - p;
              V_(TOST(i),TOST(k+1)) = V_(TOST(i),TOST(k+1)) - p * q;
            }
          }  // (s != 0)
        }  // k loop
      }  // check convergence
    }  // while (n >= low)
      
    // Backsubstitute to find vectors of upper triangular form

    if (norm == 0.0)
    {
       return;
    }
   
    for (n = nn-1; n >= 0; n--)
    {
      p = d_[TOST(n)];
      q = e_[TOST(n)];
   
      // Real vector
   
      if (q == 0)
      {
        int l = n;
        H_(TOST(n),TOST(n)) = 1.0;
        for (int i = n-1; i >= 0; i--)
        {
          w = H_(TOST(i),TOST(i)) - p;
          r = 0.0;
          for (int j = l; j <= n; j++)
          {
            r = r + H_(TOST(i),TOST(j)) * H_(TOST(j),TOST(n));
          }
          if (e_[TOST(i)] < 0.0)
          {
            z = w;
            s = r;
          }
          else
          {
            l = i;
            if (e_[TOST(i)] == 0.0)
            {
              if (w != 0.0)
              {
                H_(TOST(i),TOST(n)) = -r / w;
              }
              else
              {
                H_(TOST(i),TOST(n)) = -r / (eps * norm);
              }
   
              // Solve real equations
   
            }
            else
            {
              x = H_(TOST(i),TOST(i+1));
              y = H_(TOST(i+1),TOST(i));
              q = (d_[TOST(i)] - p) * (d_[TOST(i)] - p) + e_[TOST(i)] * e_[TOST(i)];
              t = (x * s - z * r) / q;
              H_(TOST(i),TOST(n)) = t;
              if (NumTools::abs<Real>(x) > NumTools::abs<Real>(z))
              {
                H_(TOST(i+1),TOST(n)) = (-r - w * t) / x;
              }
              else
              {
                H_(TOST(i+1),TOST(n)) = (-s - y * t) / z;
              }
            }
   
            // Overflow control
   
            t = NumTools::abs<Real>(H_(TOST(i),TOST(n)));
            if ((eps * t) * t > 1)
            {
              for (int j = i; j <= n; j++)
              {
                H_(TOST(j),TOST(n)) = H_(TOST(j),TOST(n)) / t;
              }
            }
          }
        }
   
        // Complex vector
   
      }
      else if (q < 0)
      {
        int l = n-1;

        // Last vector component imaginary so matrix is triangular
   
        if (NumTools::abs<Real>(H_(TOST(n),TOST(n-1))) > NumTools::abs<Real>(H_(TOST(n-1),TOST(n))))
        {
          H_(TOST(n-1),TOST(n-1)) = q / H_(TOST(n),TOST(n-1));
          H_(TOST(n-1),TOST(n)) = -(H_(TOST(n),TOST(n)) - p) / H_(TOST(n),TOST(n-1));
        }
        else
        {
          cdiv(0.0,-H_(TOST(n-1),TOST(n)),H_(TOST(n-1),TOST(n-1))-p,q);
          H_(TOST(n-1),TOST(n-1)) = cdivr;
          H_(TOST(n-1),TOST(n)) = cdivi;
        }
        H_(TOST(n),TOST(n-1)) = 0.0;
        H_(TOST(n),TOST(n)) = 1.0;
        for (int i = n-2; i >= 0; i--)
        {
          Real ra,sa,vr,vi;
          ra = 0.0;
          sa = 0.0;
          for (int j = l; j <= n; j++)
          {
            ra = ra + H_(TOST(i),TOST(j)) * H_(TOST(j),TOST(n-1));
            sa = sa + H_(TOST(i),TOST(j)) * H_(TOST(j),TOST(n));
          }
          w = H_(TOST(i),TOST(i)) - p;
   
          if (e_[TOST(i)] < 0.0)
          {
            z = w;
            r = ra;
            s = sa;
          }
          else
          {
            l = i;
            if (e_[TOST(i)] == 0)
            {
              cdiv(-ra,-sa,w,q);
              H_(TOST(i),TOST(n-1)) = cdivr;
              H_(TOST(i),TOST(n)) = cdivi;
            }
            else
            {
              // Solve complex equations
 
              x = H_(TOST(i),TOST(i+1));
              y = H_(TOST(i+1),TOST(i));
              vr = (d_[TOST(i)] - p) * (d_[TOST(i)] - p) + e_[TOST(i)] * e_[TOST(i)] - q * q;
              vi = (d_[TOST(i)] - p) * 2.0 * q;
              if ((vr == 0.0) && (vi == 0.0))
              {
                vr = eps * norm * (NumTools::abs<Real>(w) + NumTools::abs<Real>(q) +
                NumTools::abs<Real>(x) + NumTools::abs<Real>(y) + NumTools::abs<Real>(z));
              }
              cdiv(x*r-z*ra+q*sa,x*s-z*sa-q*ra,vr,vi);
              H_(TOST(i),TOST(n-1)) = cdivr;
              H_(TOST(i),TOST(n)) = cdivi;
              if (NumTools::abs<Real>(x) > (NumTools::abs<Real>(z) + NumTools::abs<Real>(q)))
              {
                H_(TOST(i+1),TOST(n-1)) = (-ra - w * H_(TOST(i),TOST(n-1)) + q * H_(TOST(i),TOST(n))) / x;
                H_(TOST(i+1),TOST(n)) = (-sa - w * H_(TOST(i),TOST(n)) - q * H_(TOST(i),TOST(n-1))) / x;
              }
              else
              {
                cdiv(-r-y*H_(TOST(i),TOST(n-1)),-s-y*H_(TOST(i),TOST(n)),z,q);
                H_(TOST(i+1),TOST(n-1)) = cdivr;
                H_(TOST(i+1),TOST(n)) = cdivi;
              }
            }
 
            // Overflow control
            t = std::max(NumTools::abs<Real>(H_(TOST(i),TOST(n-1))),NumTools::abs<Real>(H_(TOST(i),TOST(n))));
            if ((eps * t) * t > 1)
            {
              for (int j = i; j <= n; j++)
              {
                H_(TOST(j),TOST(n-1)) = H_(TOST(j),TOST(n-1)) / t;
                H_(TOST(j),TOST(n)) = H_(TOST(j),TOST(n)) / t;
              }
            }
          }
        }
      }
    }
   
    // Vectors of isolated roots
   
    for (int i = 0; i < nn; i++)
    {
      if (i < low || i > high)
      {
        for (int j = i; j < nn; j++)
        {
          V_(TOST(i),TOST(j)) = H_(TOST(i),TOST(j));
        }
      }
    }
   
    // Back transformation to get eigenvectors of original matrix
   
    for (int j = nn-1; j >= low; j--)
    {
      for (int i = low; i <= high; i++)
      {
        z = 0.0;
        for (int k = low; k <= std::min(j,high); k++)
        {
          z = z + V_(TOST(i),TOST(k)) * H_(TOST(k),TOST(j));
        }
        V_(TOST(i),TOST(j)) = z;
      }
    }
  }

  public:

   bool isSymmetric() const { return issymmetric_; }


    /**
     * @brief Check for symmetry, then construct the eigenvalue decomposition
     *
     * @param A    Square real (non-complex) matrix
     */
    EigenValue(const Matrix<Real>& A) :
      n_(A.getNumberOfColumns()),
      issymmetric_(true),
      d_(n_),
      e_(n_),
      V_(n_,n_),
      H_(),
      D_(n_, n_),
      ort_(),
      cdivr(), cdivi()
    {
      if (n_ > INT_MAX)
        throw Exception("EigenValue: can only be computed for matrices <= " + TextTools::toString(INT_MAX));
      for (size_t j = 0; (j < n_) && issymmetric_; j++)
      {
        for (size_t i = 0; (i < n_) && issymmetric_; i++)
        {
          issymmetric_ = (A(i,j) == A(j,i));
        }
      }

      if (issymmetric_)
      {
        for (size_t i = 0; i < n_; i++)
        {
          for (size_t j = 0; j < n_; j++)
          {
            V_(i,j) = A(i,j);
          }
        }
   
        // Tridiagonalize.
        tred2();
   
        // Diagonalize.
        tql2();
      }
      else
      {
        H_.resize(n_,n_);
        ort_.resize(n_);
         
        for (size_t j = 0; j < n_; j++)
        {
          for (size_t i = 0; i < n_; i++)
          {
            H_(i,j) = A(i,j);
          }
        }
   
        // Reduce to Hessenberg form.
        orthes();
   
        // Reduce Hessenberg to real Schur form.
        hqr2();
      }
    }


    /**
     * @brief Return the eigenvector matrix
     *
     * @return V
     */
    const RowMatrix<Real>& getV() const { return V_; }

    /**
     * @brief Return the real parts of the eigenvalues
     *
     * @return real(diag(D))
     */
    const std::vector<Real>& getRealEigenValues() const { return d_; }

    /**
     * @brief Return the imaginary parts of the eigenvalues in parameter e.
     *
     * @return e: new matrix with imaginary parts of the eigenvalues.
     */
    const std::vector<Real>& getImagEigenValues() const { return e_; }
   
    /**
     * @brief Computes the block diagonal eigenvalue matrix.
     * 
     * If the original matrix A is not symmetric, then the eigenvalue 
     * matrix D is block diagonal with the real eigenvalues in 1-by-1 
     * blocks and any complex eigenvalues,
     * a + i*b, in 2-by-2 blocks, [a, b; -b, a].  That is, if the complex
     * eigenvalues look like
     * <pre>
     *
     *       u + iv     .        .          .      .    .
     *         .      u - iv     .          .      .    .
     *         .        .      a + ib       .      .    .
     *         .        .        .        a - ib   .    .
     *         .        .        .          .      x    .
     *         .        .        .          .      .    y
     * </pre>
     *     then D looks like
     * <pre>
     *
     *         u        v        .          .      .    .
     *        -v        u        .          .      .    . 
     *         .        .        a          b      .    .
     *         .        .       -b          a      .    .
     *         .        .        .          .      x    .
     *         .        .        .          .      .    y
     * </pre>
     * This keeps V a real matrix in both symmetric and non-symmetric
     * cases, and A*V = V*D.
     *
     * @return D: upon return, the matrix is filled with the block diagonal 
     * eigenvalue matrix.
     */
    const RowMatrix<Real>& getD() const
    {
      for (size_t i = 0; i < n_; i++)
      {
        for (size_t j = 0; j < n_; j++)
        {
          D_(i,j) = 0.0;
        }
        D_(i,i) = d_[i];
        if (e_[i] > 0)
        {
          D_(i,i+1) = e_[i];
        }
        else if (e_[i] < 0)
        {
          D_(i,i-1) = e_[i];
        }
      }
      return D_;
    }
};

} //end of namespace bpp.

#endif //_EIGENVALUE_H_

