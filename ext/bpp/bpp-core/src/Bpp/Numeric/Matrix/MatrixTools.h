//
// File: MatrixTools.h
// Created by: Julien Dutheil
// Created on: Mon Jan 19 16:42:25 2004
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

#ifndef _MATRIXTOOLS_H_
#define _MATRIXTOOLS_H_

#include "../VectorTools.h"
#include "Matrix.h"
#include "LUDecomposition.h"
#include "EigenValue.h"
#include "../../Io/OutputStream.h"

#include <cstdio>
#include <iostream>

namespace bpp
{
/**
 * @brief Functions dealing with matrices.
 */
  class MatrixTools
  {
  public:
    MatrixTools() {}
    ~MatrixTools() {}

  public:
    /**
     * @brief Copy operation. This function supplies the lack of inheritence of the assigment operator :D .
     *
     * @param A [in] Original matrix.
     * @param O [out] A copy of the given matrix.
     */
    template<class MatrixA, class MatrixO>
    static void copy(const MatrixA& A, MatrixO& O)
    {
      O.resize(A.getNumberOfRows(), A.getNumberOfColumns());
      for (size_t i = 0; i < A.getNumberOfRows(); i++)
      {
        for (size_t j = 0; j < A.getNumberOfColumns(); j++)
        {
          O(i, j) = A(i, j);
        }
      }
    }

    /**
     * @brief Get a identity matrix of a given size.
     *
     * @param n the size of the matrix.
     * @param O [out] A identity matrix of size n.
     */
    template<class Matrix>
    static void getId(size_t n, Matrix& O)
    {
      O.resize(n, n);
      for (size_t i = 0; i < n; i++)
      {
        for (size_t j = 0; j < n; j++) {
          O(i, j) = (i == j) ? 1 : 0;
        }
      }
    }

    /**
     * @param D [in] A vector of diagonal elements.
     * @param O [out] A diagonal matrix with diagonal elements taken from a vector.
     */
    template<class Scalar>
    static void diag(const std::vector<Scalar>& D, Matrix<Scalar>& O)
    {
      size_t n = D.size();
      O.resize(n, n);
      for (size_t i = 0; i < n; i++)
      {
        for (size_t j = 0; j < n; j++) { O(i, j) = (i == j) ? D[i] : 0;}
      }
    }

    /**
     * @param x [in] A scalar
     * @param n [in] the dimension of the output matrix
     * @param O [out] A diagonal matrix with diagonal elements equal to x
     */
    template<class Scalar>
    static void diag(const Scalar x, size_t n, Matrix<Scalar>& O)
    {
      O.resize(n, n);
      for (size_t i = 0; i < n; i++)
      {
        for (size_t j = 0; j < n; j++) { O(i, j) = (i == j) ? x : 0;}
      }
    }

    /**
     * @param M [in] The matrix.
     * @param O [out] The diagonal elements of a square matrix as a vector.
     * @throw DimensionException If M is not a square matrix.
     */
    template<class Scalar>
    static void diag(const Matrix<Scalar>& M, std::vector<Scalar>& O)
    noexcept(false)
    {
      size_t nc = M.getNumberOfColumns();
      size_t nr = M.getNumberOfRows();
      if (nc != nr) throw DimensionException("MatrixTools::diag(). M must be a square matrix.", nr, nc);
      O.resize(nc);
      for (size_t i = 0; i < nc; i++) { O[i] = M(i, i);}
    }

    /**
     * @brief Set all elements in M to value x.
     * @param M A matrix.
     * @param x The value to use.
     */
    template<class Matrix, class Scalar>
    static void fill(Matrix& M, Scalar x)
    {
      for (size_t i = 0; i < M.getNumberOfRows(); i++)
      {
        for (size_t j = 0; j < M.getNumberOfColumns(); j++)
        {
          M(i, j) = x;
        }
      }
    }

    /**
     * @brief Set all diagonal elements in M to value x.
     * @param M A matrix.
     * @param x The value to use.
     */
    template<class Matrix, class Scalar>
    static void fillDiag(Matrix& M, Scalar x)
    {
      for (size_t i = 0; i < M.getNumberOfRows(); i++)
        M(i, i) = x;
    }

    /**
     * @brief Multiply all elements of a matrix by a given value, and add a constant.
     *
     * Performs \f$\forall i \forall j m_{i,j} = a.m_{i,j}+b\f$.
     *
     * @param A A matrix.
     * @param a Multiplicator.
     * @param b Constant.
     */
    template<class Matrix, class Scalar>
    static void scale(Matrix& A, Scalar a, Scalar b = 0)
    {
      for (size_t i = 0; i < A.getNumberOfRows(); i++)
      {
        for (size_t j = 0; j < A.getNumberOfColumns(); j++)
        {
          A(i, j) = a * A(i, j) + b;
        }
      }
    }

    /**
     * @param A [in] First matrix.
     * @param B [in] Second matrix.
     * @param O [out] The dot product of two matrices.
     */
    template<class Scalar>
    static void mult(const Matrix<Scalar>& A, const Matrix<Scalar>& B, Matrix<Scalar>& O)
    noexcept(false)
    {
      size_t ncA = A.getNumberOfColumns();
      size_t nrA = A.getNumberOfRows();
      size_t nrB = B.getNumberOfRows();
      size_t ncB = B.getNumberOfColumns();
      if (ncA != nrB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrB, ncA);
      O.resize(nrA, ncB);
      for (size_t i = 0; i < nrA; i++)
      {
        for (size_t j = 0; j < ncB; j++)
        {
          O(i, j) = 0;
          for (size_t k = 0; k < ncA; k++)
          {
            O(i, j) += A(i, k) * B(k, j);
          }
        }
      }
    }

    /**
     * @brief Compute A . D . B where D is a diagonal matrix in O(n^3).
     *
     * Since D is a diagonal matrix, this function is more efficient than doing
     * mult(mult(A, diag(D)), B), which involves two 0(n^3) operations.
     *
     * @param A [in] The first matrix.
     * @param D [in] The diagonal matrix (only diagonal elements in a vector)
     * @param B [in] The second matrix.
     * @param O [out] The result matrix.
     * @throw DimensionException If matrices have not the appropriate size.
     */
    template<class Scalar>
    static void mult(const Matrix<Scalar>& A, const std::vector<Scalar>& D, const Matrix<Scalar>& B, Matrix<Scalar>& O)
    noexcept(false)
    {
      size_t ncA = A.getNumberOfColumns();
      size_t nrA = A.getNumberOfRows();
      size_t nrB = B.getNumberOfRows();
      size_t ncB = B.getNumberOfColumns();
      if (ncA != nrB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrB, ncA);
      if (ncA != D.size()) throw DimensionException("MatrixTools::mult(). Vector size is not equal to matrix size.", D.size(), ncA);
      O.resize(nrA, ncB);
      for (size_t i = 0; i < nrA; i++)
      {
        for (size_t j = 0; j < ncB; j++)
        {
          O(i, j) = 0;
          for (size_t k = 0; k < ncA; k++)
          {
            O(i, j) += A(i, k) * B(k, j) * D[k];
          }
        }
      }
    }

    /**
     * @brief Compute A . (U+D+L) . B where D is a diagonal matrix, U
     * (resp. L) is a matrix in which the only non-zero terms are on the
     * diagonal that is over (resp. under) the main diagonal, in O(n^3).
     *
     * Since D is a diagonal matrix, this function is more efficient than doing
     * mult(mult(A, diag(D)), B), which involves two 0(n^3) operations.
     *
     * @param A [in] The first matrix.
     * @param D [in] The diagonal matrix (only diagonal elements in a vector)
     * @param U [in] The upper diagonal matrix (only upper diagonal elements in a vector)
     * @param L [in] The lower diagonal matrix (only lower diagonal elements in a vector)
     * @param B [in] The second matrix.
     * @param O [out] The result matrix.
     * @throw DimensionException If matrices have not the appropriate size.
     */
    template<class Scalar>
    static void mult(const Matrix<Scalar>& A, const std::vector<Scalar>& D, const std::vector<Scalar>& U, const std::vector<Scalar>& L, const Matrix<Scalar>& B, Matrix<Scalar>& O)
    noexcept(false)
    {
      size_t ncA = A.getNumberOfColumns();
      size_t nrA = A.getNumberOfRows();
      size_t nrB = B.getNumberOfRows();
      size_t ncB = B.getNumberOfColumns();
      if (ncA != nrB) throw DimensionException("MatrixTools::mult(). nrows B != ncols A.", nrB, ncA);
      if (ncA != D.size()) throw DimensionException("MatrixTools::mult(). Vector size is not equal to matrix size.", D.size(), ncA);
      if (ncA != U.size()+1) throw DimensionException("MatrixTools::mult(). Vector size is not equal to matrix size-1.", U.size(), ncA);
      if (ncA != L.size()+1) throw DimensionException("MatrixTools::mult(). Vector size is not equal to matrix size-1.", L.size(), ncA);
      O.resize(nrA, ncB);
      for (size_t i = 0; i < nrA; i++)
      {
        for (size_t j = 0; j < ncB; j++)
        {
          O(i, j) = A(i, 0) * D[0] * B(0, j);
          if (nrB>1)
            O(i, j) += A(i,0) * U[0] * B(1,j);
          for (size_t k = 1; k < ncA-1; k++)
          {
            O(i, j) += A(i, k) * (L[k-1] * B(k-1, j) + D[k] * B(k, j) + U[k] * B(k+1,j));
          }
          if (ncA>=2)
            O(i, j) += A(i, ncA-1) * L[ncA-2] * B(ncA-2, j);
          O(i,j) += A(i, ncA-1) * D[ncA-1] * B(ncA-1, j);
        }
      }
    }

    /**
     * @brief Add matrix B to matrix A.
     *
     * @param A [in, out] Matrix A
     * @param B [in] Matrix B
     * @throw DimensionException If A and B have note the same size.
     */
    
    template<class MatrixA, class MatrixB>
    static void add(MatrixA& A, const MatrixB& B) noexcept(false)
    {
      size_t ncA = A.getNumberOfColumns();
      size_t nrA = A.getNumberOfRows();
      size_t nrB = B.getNumberOfRows();
      size_t ncB = B.getNumberOfColumns();
      if (ncA > ncB) throw DimensionException("MatrixTools::operator+=(). A and B must have the same number of colums.", ncB, ncA);
      if (nrA > nrB) throw DimensionException("MatrixTools::operator+=(). A and B must have the same number of rows.", nrB, nrA);
      
      
      for (size_t i = 0; i < nrA; i++)
      {
        for (size_t j = 0; j < ncA; j++)
        {
          A(i, j) += B(i, j);
        }
      }
    }

    /**
     * @brief Add matrix x.B to matrix A.
     *
     * @param A [in,out] Matrix A
     * @param x [in] Scalar x
     * @param B [in] Matrix B
     * @throw DimensionException If A and B have note the same size.
     */
    template<class MatrixA, class MatrixB, class Scalar>
    static void add(MatrixA& A, Scalar& x, const MatrixB& B) noexcept(false)
    {
      size_t ncA = A.getNumberOfColumns();
      size_t nrA = A.getNumberOfRows();
      size_t nrB = B.getNumberOfRows();
      size_t ncB = B.getNumberOfColumns();
      if (ncA != ncB) throw DimensionException("MatrixTools::operator+(). A and B must have the same number of colums.", ncB, ncA);
      if (nrA != nrB) throw DimensionException("MatrixTools::operator+(). A and B must have the same number of rows.", nrB, nrA);
      
      for (size_t i = 0; i < nrA; i++)
      {
        for (size_t j = 0; j < ncA; j++)
        {
          A(i, j) += x*B(i, j);
        }
      }
    }

    /**
     * @brief Compute the power of a given matrix.
     *
     * @param A [in] The matrix.
     * @param p The number of multiplications.
     * @param O [out]\f$  A^p \f$ computed recursively:
     *               \f$ A^{2n} = (A^n)^2 \f$
     *               \f$ A^{2n+1} = A*(A^n)^2 \f$   
     * If p = 0, sends the identity matrix.
     * @throw DimensionException If m is not a square matrix.
     */
    template<class Matrix>
    static void pow(const Matrix& A, size_t p, Matrix& O) noexcept(false)
    {
      size_t n = A.getNumberOfRows();
      if (n != A.getNumberOfColumns()) throw DimensionException("MatrixTools::pow(). nrows != ncols.", A.getNumberOfColumns(), A.getNumberOfRows());
      switch(p){
      case 0:
        getId<Matrix>(n, O);
        break;
      case 1:
        copy(A,O);
        break;
      case 2:
        mult(A,A,O);
        break;
      default:
        Matrix tmp;
        if (p%2){
          pow(A,p/2,tmp);
          pow(tmp,2,O);
        }
        else{
          pow(A,(p-1)/2,tmp);
          pow(tmp,2,O);
          mult(A,O,tmp);
          copy(tmp,O);
        }
      }
    }

    /**
     * @brief Compute the power of a given matrix, using eigen value decomposition.
     *
     * @param A [in] The matrix.
     * @param p The power of the matrix.
     * @param O [out]\f$\prod_{i=1}^p m\f$.
     * If p = 0, sends the identity matrix.
     * @throw DimensionException If m is not a square matrix.
     */
    template<class Scalar>
    static void pow(const Matrix<Scalar>& A, double p, Matrix<Scalar>& O) noexcept(false)
    {
      size_t n = A.getNumberOfRows();
      if (n != A.getNumberOfColumns()) throw DimensionException("MatrixTools::pow(). nrows != ncols.", A.getNumberOfColumns(), A.getNumberOfRows());
      EigenValue<Scalar> eigen(A);
      RowMatrix<Scalar> rightEV, leftEV;
      rightEV = eigen.getV();
      inv(rightEV, leftEV);
      mult(rightEV, VectorTools::pow(eigen.getRealEigenValues(), p), leftEV, O);
    }

    /**
     * @brief Perform matrix exponentiation using diagonalization.
     *
     * @warning This method currently relies only on diagonalization, so it won't work if your matrix is not diagonalizable.
     * The function may be extended later to deal with other cases.
     *
     * @param A [in] The matrix.
     * @param O [out]\f$\prod_{i=1}^p m\f$.
     * @throw DimensionException If m is not a square matrix.
     */
    template<class Scalar>
    static void exp(const Matrix<Scalar>& A, Matrix<Scalar>& O) noexcept(false)
    {
      size_t n = A.getNumberOfRows();
      if (n != A.getNumberOfColumns()) throw DimensionException("MatrixTools::exp(). nrows != ncols.", A.getNumberOfColumns(), A.getNumberOfRows());
      EigenValue<Scalar> eigen(A);
      RowMatrix<Scalar> rightEV, leftEV;
      rightEV = eigen.getV();
      inv(rightEV, leftEV);
      mult(rightEV, VectorTools::exp(eigen.getRealEigenValues()), leftEV, O);
    }

    /**
     * @brief Compute a vector of the first powers of a given matrix.
     *
     * @param A [in] The matrix.
     * @param p The number of powers.
     * @param vO [out] the vector of the powers (from 0 to p)
     *
     * @throw DimensionException If m is not a square matrix.
     */
  
    template<class Matrix, class Scalar>
    static void Taylor(const Matrix& A, size_t p, std::vector< RowMatrix<Scalar> > & vO)
    noexcept(false)
    {
      size_t n = A.getNumberOfRows();
      if (n != A.getNumberOfColumns())
        throw DimensionException("MatrixTools::pow(). nrows != ncols.", A.getNumberOfColumns(), A.getNumberOfRows());
      vO.resize(p+1);
      getId<Matrix>(n, vO[0]);
      copy(A,vO[1]);
    
      for (size_t i = 1; i < p; i++)
      {
        mult(vO[i], A, vO[i+1]);
      }
    }

    /**
     * @return The position of the maximum value in the matrix.
     * @param m The matrix.
     */
    template<class Matrix>
    static std::vector<size_t> whichMax(const Matrix& m)
    {
      size_t nrows = m.getNumberOfRows();
      size_t ncols = m.getNumberOfColumns();
      std::vector<size_t> pos(2);
      size_t imax = 0;
      size_t jmax = 0;
      double currentMax = log(0.);
      for (size_t i = 0; i < nrows; i++)
      {
        for (size_t j = 0; j < ncols; j++)
        {
          double currentValue = m(i, j);
          // cout << currentValue << "\t" << (currentValue > currentMax) << endl;
          if (currentValue > currentMax)
          {
            imax = i;
            jmax = j;
            currentMax = currentValue;
          }
        }
      }
      pos[0] = imax;
      pos[1] = jmax;
      return pos;
    }

    /**
     * @return The position of the minimum value in the matrix.
     * @param m The matrix.
     */
    template<class Matrix>
    static std::vector<size_t> whichMin(const Matrix& m)
    {
      size_t nrows = m.getNumberOfRows();
      size_t ncols = m.getNumberOfColumns();
      std::vector<size_t> pos(2);
      size_t imin = 0;
      size_t jmin = 0;
      double currentMin = -log(0.);
      for (size_t i = 0; i < nrows; i++)
      {
        for (size_t j = 0; j < ncols; j++)
        {
          double currentValue = m(i, j);
          if (currentValue < currentMin)
          {
            imin = i;
            jmin = j;
            currentMin = currentValue;
          }
        }
      }
      pos[0] = imin;
      pos[1] = jmin;
      return pos;
    }

    /**
     * @return The maximum value in the matrix.
     * @param m The matrix.
     */
    template<class Real>
    static Real max(const Matrix<Real>& m)
    {
      size_t nrows = m.getNumberOfRows();
      size_t ncols = m.getNumberOfColumns();
      Real currentMax = log(0.);
      for (size_t i = 0; i < nrows; i++)
      {
        for (size_t j = 0; j < ncols; j++)
        {
          Real currentValue = m(i, j);
          // cout << currentValue << "\t" << (currentValue > currentMax) << endl;
          if (currentValue > currentMax)
          {
            currentMax = currentValue;
          }
        }
      }
      return currentMax;
    }


    /**
     * @return The minimum value in the matrix.
     * @param m The matrix.
     */
    template<class Real>
    static Real min(const Matrix<Real>& m)
    {
      size_t nrows = m.getNumberOfRows();
      size_t ncols = m.getNumberOfColumns();
      Real currentMin = -log(0.);
      for (size_t i = 0; i < nrows; i++)
      {
        for (size_t j = 0; j < ncols; j++)
        {
          Real currentValue = m(i, j);
          if (currentValue < currentMin)
          {
            currentMin = currentValue;
          }
        }
      }
      return currentMin;
    }

    /**
     * @brief Print a matrix to a stream.
     *
     * @param m The matrix to print.
     * @param out The stream to use.
     */
    template<class Matrix>
    static void print(const Matrix& m, std::ostream& out = std::cout)
    {
      out << m.getNumberOfRows() << "x" << m.getNumberOfColumns() << std::endl;
      out << "[" << std::endl;
      for (size_t i = 0; i < m.getNumberOfRows(); i++)
      {
        out << "[";
        for (size_t j = 0; j < m.getNumberOfColumns() - 1; j++)
        {
          out << m(i, j) << ", ";
        }
        if (m.getNumberOfColumns() > 0) out << m(i, m.getNumberOfColumns() - 1) << "]" << std::endl;
      }
      out << "]" << std::endl;
    }

    /**
     * @brief Print a matrix to a stream.
     *
     * @param m The matrix to print.
     * @param out The stream to use.
     * @param pIn left delimiter (default: "(")
     * @param pOut right delimiter (default: ")")
     */
    template<class Matrix>
    static void print(const Matrix& m, bpp::OutputStream& out, char pIn = '(', char pOut = ')')
    {
      out << pIn;
    
      for (size_t i = 0; i < m.getNumberOfRows(); i++)
      {
        if (i!=0)
          out << ",";
      
        out << pIn;
        for (size_t j = 0; j < m.getNumberOfColumns() - 1; j++)
        {
          out << m(i, j) << ", ";
        }
        if (m.getNumberOfColumns() > 0) out << m(i, m.getNumberOfColumns() - 1) << pOut;
      }
      out << pOut;
    }

    /**
     * @brief Print a matrix to a stream, so that it is read by R.
     *
     * @param m The matrix to print.
     * @param variableName The name of the R variable handeling the matrix
     * @param out The stream to use.
     */
    template<class Matrix>
    static void printForR(const Matrix& m, const std::string& variableName = "x", std::ostream& out = std::cout)
    {
      out.precision(12);
      out << variableName << "<-matrix(c(";
      for (size_t i = 0; i < m.getNumberOfRows(); i++)
      {
        for (size_t j = 0; j < m.getNumberOfColumns(); j++)
        {
          if (i > 0 || j > 0)
            out << ", ";
          out << m(i, j);
        }
      }
      out << "), nrow=" << m.getNumberOfRows() << ", byrow=TRUE)" << std::endl;
    }


    /**
     * @brief Print a vector to a stream.
     *
     * @param v The vector to print.
     * @param out The stream to use.
     */
    template<class Real>
    static void print(const std::vector<Real>& v, std::ostream& out = std::cout)
    {
      out << v.size() << std::endl;
      out << "[";
      for (size_t i = 0; i < v.size() - 1; i++)
      {
        out << v[i] << ", ";
      }
      if (v.size() > 0) out << v[v.size() - 1];
      out << "]" << std::endl;
    }

    /**
     * @return True if the matrix is a square matrix.
     * @param A A matrix.
     */
    template<class Matrix>
    static bool isSquare(const Matrix& A) { return A.getNumberOfRows() == A.getNumberOfColumns(); }

    /**
     * @param A [in] The matrix to inverse.
     * @param O [out] The inverse matrix of A.
     * @return x the minimum absolute value of the diagonal of the LU decomposition
     * @throw DimensionException If A is not a square matrix.
     */
    template<class Scalar>
    static Scalar inv(const Matrix<Scalar>& A, Matrix<Scalar>& O)
    noexcept(false)
    {
      if (!isSquare(A)) throw DimensionException("MatrixTools::inv(). Matrix A is not a square matrix.", A.getNumberOfRows(), A.getNumberOfColumns());
      LUDecomposition<Scalar> lu(A);
      RowMatrix<Scalar> I;
      getId(A.getNumberOfRows(), I);
      return lu.solve(I, O);
    }

    /**
     * @brief Get determinant of a square matrix.
     *
     * This implementation is in @f$o(n^3)@f$ and uses the LU decomposition method.
     *
     * @param A [in] The input matrix.
     * @return The determinant of A.
     * @throw DimensionException If A is not a square matrix.
     */
    template<class Scalar>
    static double det(const Matrix<Scalar>& A) noexcept(false)
    {
      if (!isSquare(A)) throw DimensionException("MatrixTools::det(). Matrix A is not a square matrix.", A.getNumberOfRows(), A.getNumberOfColumns());
      LUDecomposition<Scalar> lu(A);
      return lu.det();
    }

    /**
     * @param A [in] The matrix to transpose.
     * @param O [out] The transposition of A.
     */
    template<class MatrixA, class MatrixO>
    static void transpose(const MatrixA& A, MatrixO& O)
    {
      O.resize(A.getNumberOfColumns(), A.getNumberOfRows());
      for (size_t i = 0; i < A.getNumberOfColumns(); i++)
      {
        for (size_t j = 0; j < A.getNumberOfRows(); j++)
        {
          O(i, j) = A(j, i);
        }
      }
    }

    /**
     * @brief Compute the variance-covariance matrix of an input matrix.
     *
     * The input matrix represent a n-sample of a random vector of dimension r.
     * It is assumed to have r rows and n columns.
     * The variance matrix is then computed as @f[ V = A\cdot A^T - \mu\cdot\mu^T@f],
     * where @f$\mu@f$ is the mean vector of the sample.
     * the output matrix is a square matrix of size r.
     *
     * @param A [in] The intput matrix.
     * @param O [out] The resulting variance covariance matrix.
     */
    template<class Scalar>
    static void covar(const Matrix<Scalar>& A, Matrix<Scalar>& O)
    {
      size_t r = A.getNumberOfRows();
      size_t n = A.getNumberOfColumns();
      O.resize(r, r);
      RowMatrix<Scalar> tA;
      transpose(A, tA);
      mult(A, tA, O);
      scale(O, 1. / static_cast<double>(n));
      RowMatrix<Scalar> mean(r, 1);
      for (size_t i = 0; i < r; i++)
      {
        for (size_t j = 0; j < n; j++)
        {
          mean(i, 0) += A(i, j);
        }
        mean(i, 0) /= static_cast<double>(n);
      }
      RowMatrix<Scalar> tMean;
      transpose(mean, tMean);
      RowMatrix<Scalar> meanMat;
      mult(mean, tMean, meanMat);
      scale(meanMat, -1.);
      add(O, meanMat);
    }

    /**
     * @brief Compute the Kronecker product of two row matrices.
     *
     * @param A [in] The first row matrix.
     * @param B [in] The second row matrix.
     * @param O [out] The product \f$A \otimes B\f$.
     * @param check [optional] if resize of 0 (default: true)
     */
    template<class Scalar>
    static void kroneckerMult(const Matrix<Scalar>& A, const Matrix<Scalar>& B, Matrix<Scalar>& O, bool check = true)
    {
      size_t ncA = A.getNumberOfColumns();
      size_t nrA = A.getNumberOfRows();
      size_t nrB = B.getNumberOfRows();
      size_t ncB = B.getNumberOfColumns();

      if (check)
        O.resize(nrA * nrB, ncA * ncB);
      
      for (size_t ia = 0; ia < nrA; ia++)
      {
        for (size_t ja = 0; ja < ncA; ja++)
        {
          Scalar aij = A(ia, ja);
          for (size_t ib = 0; ib < nrB; ib++)
          {
            for (size_t jb = 0; jb < ncB; jb++)
            {
              O(ia * nrB + ib, ja * ncB + jb) = aij * B(ib, jb);
            }
          }
        }
      }
    }

    /**
     * @brief Compute the Kronecker product of one Matrice with a
     * diagonal matrix, which main term and dimension are given
     *
     * @param A [in] The first row matrix.
     * @param dim [in] The dimension of the diagonal matrix.
     * @param v [in] The diagonal value of the diagonal matrix.
     * @param O [out] The product \f$A \otimes B\f$.
     * @param check [optional] if resize of 0 (default: true)
     */
    template<class Scalar>
    static void kroneckerMult(const Matrix<Scalar>& A, size_t dim, const Scalar& v, Matrix<Scalar>& O, bool check = true)
    {
      size_t ncA = A.getNumberOfColumns();
      size_t nrA = A.getNumberOfRows();

      if (check)
        O.resize(nrA * dim, ncA * dim);
      
      for (size_t ia = 0; ia < nrA; ia++)
      {
        for (size_t ja = 0; ja < ncA; ja++)
        {
          Scalar aij = A(ia, ja);
          for (size_t ib = 0; ib < dim; ib++)
          {
            for (size_t jb = 0; jb < dim; jb++)
            {
              O(ia * dim + ib, ja * dim + jb) = aij * ((ib==jb)?v:0);
            }
          }
        }
      }
    }

    /**
     * @brief Compute the Kronecker product of two row matrices in
     * which the diagonal element is changed
     *
     * @param A [in] The first row matrix.
     * @param B [in] The second row matrix.
     * @param dA [in] The replaced diagonal element of A
     * @param dB [in] The replaced diagonal element of B
     * @param O [out] The product \f$A \otimes B\f$.
     * @param check [optional] if resize of 0 (default: true)
     */
    template<class Scalar>
    static void kroneckerMult(const Matrix<Scalar>& A, const Matrix<Scalar>& B, const Scalar& dA, const Scalar& dB, Matrix<Scalar>& O, bool check= true)
    {
      size_t ncA = A.getNumberOfColumns();
      size_t nrA = A.getNumberOfRows();
      size_t nrB = B.getNumberOfRows();
      size_t ncB = B.getNumberOfColumns();

      if (check)
        O.resize(nrA * nrB, ncA * ncB);
      
      for (size_t ia = 0; ia < nrA; ia++)
      {
        for (size_t ja = 0; ja < ncA; ja++)
        {
          const Scalar& aij = (ia==ja)?dA:A(ia, ja);
          
          for (size_t ib = 0; ib < nrB; ib++)
          {
            for (size_t jb = 0; jb < ncB; jb++)
            {
              O(ia * nrB + ib, ja * ncB + jb) = aij * ((ib==jb)?dB:B(ib, jb));
            }
          }
        }
      }
    }

    /**
     * @brief Compute the Hadamard product of two row matrices with same dimensions.
     *
     * @param A [in] The first row matrix.
     * @param B [in] The second row matrix.
     * @param O [out] The Hadamard product.
     */
    template<class Scalar>
    static void hadamardMult(const Matrix<Scalar>& A, const Matrix<Scalar>& B, Matrix<Scalar>& O)
    {
      size_t ncA = A.getNumberOfColumns();
      size_t nrA = A.getNumberOfRows();
      size_t nrB = B.getNumberOfRows();
      size_t ncB = B.getNumberOfColumns();
      if (nrA != nrB) throw DimensionException("MatrixTools::hadamardMult(). nrows A != nrows B.", nrA, nrB);
      if (ncA != ncB) throw DimensionException("MatrixTools::hadamardMult(). ncols A != ncols B.", ncA, ncB);
      O.resize(nrA, ncA);
      for (size_t i = 0; i < nrA; i++)
      {
        for (size_t j = 0; j < ncA; j++)
        {
          O(i, j) = A(i, j) * B(i, j);
        }
      }
    }

    /**
     * @brief Compute the "Hadamard" product of a row matrix and a vector containing weights, according to rows or columns.
     *
     * @param A [in] The row matrix.
     * @param B [in] The vector of row or column weights.
     * @param O [out] The 'Hadamard' product.
     * @param row Boolean. If row is set to 'true', the vector contains weights for rows. Otherwise the vector contains weights for columns.
     */
    template<class Scalar>
    static void hadamardMult(const Matrix<Scalar>& A, const std::vector<Scalar>& B, Matrix<Scalar>& O, bool row = true)
    {
      size_t ncA = A.getNumberOfColumns();
      size_t nrA = A.getNumberOfRows();
      size_t sB = B.size();
      if (row == true && nrA != sB) throw DimensionException("MatrixTools::hadamardMult(). nrows A != size of B.", nrA, sB);
      if (row == false && ncA != sB) throw DimensionException("MatrixTools::hadamardMult(). ncols A != size of B.", ncA, sB);
      O.resize(nrA, ncA);
      if (row)
      {
        for (size_t i = 0; i < nrA; i++)
        {
          for (size_t j = 0; j < ncA; j++)
          {
            O(i, j) = A(i, j) * B[i];
          }
        }
      }
      else
      {
        for (size_t i = 0; i < nrA; i++)
        {
          for (size_t j = 0; j < ncA; j++)
          {
            O(i, j) = A(i, j) * B[j];
          }
        }
      }
    }

    /**
     * @brief Compute the direct sum of two row matrices.
     *
     * @param A [in] The first row matrix.
     * @param B [in] The second row matrix.
     * @param O [out] The sum \f$A \oplus B\f$.
     */
    template<class Scalar>
    static void directSum(const Matrix<Scalar>& A, const Matrix<Scalar>& B, Matrix<Scalar>& O)
    {
      size_t ncA = A.getNumberOfColumns();
      size_t nrA = A.getNumberOfRows();
      size_t nrB = B.getNumberOfRows();
      size_t ncB = B.getNumberOfColumns();
      O.resize(nrA + nrB, ncA + ncB);
      
      for (size_t ia = 0; ia < nrA; ia++)
      {
        for (size_t ja = 0; ja < ncA; ja++)
        {
          O(ia, ja) = A(ia, ja);
        }
      }

      for (size_t ia = 0; ia < nrA; ia++)
      {
        for (size_t jb = 0; jb < ncB; jb++)
        {
          O(ia, ncA + jb) = 0;
        }
      }

      for (size_t ib = 0; ib < nrB; ib++)
      {
        for (size_t ja = 0; ja < ncA; ja++)
        {
          O(nrA + ib, ja) = 0;
        }
      }

      for (size_t ib = 0; ib < nrB; ib++)
      {
        for (size_t jb = 0; jb < nrB; jb++)
        {
          O(nrA + ib, ncA + jb) = B(ib, jb);
        }
      }
    }

    /**
     * @brief Compute the direct sum of n row matrices.
     *
     * @param vA [in] A vector of row matrices of any size.
     * @param O [out] The sum \f$\bigoplus_i A_i\f$.
     */
    template<class Scalar>
    static void directSum(const std::vector< Matrix<Scalar>*>& vA, Matrix<Scalar>& O)
    {
      size_t nr = 0;
      size_t nc = 0;
      for (size_t k = 0; k < vA.size(); k++)
      {
        nr += vA[k]->getNumberOfRows();
        nc += vA[k]->getNumberOfColumns();
      }
      O.resize(nr, nc);
      for (size_t k=0; k<nr; k++)
        for (size_t k2=0; k2<nc; k2++)
          O(k,k2)=0;
      
      size_t rk = 0; // Row counter
      size_t ck = 0; // Col counter
      for (size_t k = 0; k < vA.size(); k++)
      {
        const Matrix<Scalar>* Ak = vA[k];
        for (size_t i = 0; i < Ak->getNumberOfRows(); i++)
        {
          for (size_t j = 0; j < Ak->getNumberOfColumns(); j++)
          {
            O(rk + i, ck + j) = (*Ak)(i, j);
          }
        }
        rk += Ak->getNumberOfRows();
        ck += Ak->getNumberOfColumns();
      }
    }

    /**
     * @brief Convert to a vector of vector.
     *
     * @param M [in] A matrix object.
     * @param vO [out] The output vector of vector (will be resized accordingly).
     */
    template<class Scalar>
    static void toVVdouble(const Matrix<Scalar>& M, std::vector< std::vector<Scalar> >& vO)
    {
      size_t n = M.getNumberOfRows();
      size_t m = M.getNumberOfColumns();
      vO.resize(n);
      for (size_t i = 0; i < n; i++)
      {
        vO[i].resize(m);
        for (size_t j = 0; j < m; j++)
        {
          vO[i][j] = M(i, j);
        }
      }
    }

    /**
     * @brief Sum all elements in M.
     * @param M A matrix.
     * @return The sum of all elements.
     */
    template<class Scalar>
    static Scalar sumElements(const Matrix<Scalar>& M)
    {
      Scalar sum = 0;
      for (size_t i = 0; i < M.getNumberOfRows(); i++)
      {
        for (size_t j = 0; j < M.getNumberOfColumns(); j++)
        {
          sum += M(i, j);
        }
      }
      return sum;
    }


  
    /**
     * @brief Linear Assignment Problem
     *
     * The algorithm coded here is described in 
     * * A Shortest Augmenting Path Algorithm for Dense and Sparse Linear Assignment Problems, Computing 38, 325-340, 1987
     * by R. Jonker and A. Volgenant, University of Amsterdam.
     *
     * @param assignCost [input/output] Cost matrix
     * @param rowSol     [output] Column assigned to row in solution
     * @param colSol     [output] Row assigned to column in solution
     * @param u          [output] Dual variables, row reduction numbers
     * @param v          [output] Dual variables, column reduction numbers
     * @return The optimal cost.
     */
    template<class Scalar>
    static Scalar lap(Matrix<Scalar>& assignCost,
                      std::vector<int> &rowSol, 
                      std::vector<int> &colSol, 
                      std::vector<Scalar> &u, 
                      std::vector<Scalar> &v)
    noexcept(false)
    {
      size_t dim = assignCost.getNumberOfRows();
      if (assignCost.getNumberOfColumns() != dim)
        throw Exception("MatrixTools::lap. Cost matrix should be scare.");
  
      bool unassignedFound;
      size_t i, iMin;
      size_t numFree = 0, previousNumFree, f, k, freeRow;
      int i0;
      std::vector<size_t> free(dim); // list of unassigned rows.
      std::vector<size_t> pred(dim); // row-predecessor of column in augmenting/alternating path.
      size_t j, j1, j2, endOfPath, last, low, up;
      std::vector<size_t> colList(dim); // list of columns to be scanned in various ways.
      std::vector<short int> matches(dim, 0); // counts how many times a row could be assigned.
      Scalar min;
      Scalar h;
      size_t uMin, uSubMin;
      Scalar v2;
      std::vector<Scalar> d(dim); // 'cost-distance' in augmenting path calculation.

      // Column reduction
      for (j = dim; j > 0; j--)    // reverse order gives better results.
      {
        // find minimum cost over rows.
        min = assignCost(0, j - 1); 
        iMin = 0;
        for (i = 1; i < dim; ++i) { 
          if (assignCost(i, j - 1) < min) 
          { 
            min = assignCost(i, j - 1); 
            iMin = i;
          }
        }
        v[j - 1] = min; 

        if (++matches[iMin] == 1) 
        { 
          // init assignment if minimum row assigned for first time.
          rowSol[iMin] = static_cast<int>(j - 1); 
          colSol[j - 1] = static_cast<int>(iMin); 
        }
        else
          colSol[j - 1] = -1;        // row already assigned, column not assigned.
      }

      // Reduction tranfer
      for (i = 0; i < dim; i++) { 
        if (matches[i] == 0)     // fill list of unassigned 'free' rows.
          free[numFree++] = i;
        else {
          if (matches[i] == 1)   // transfer reduction from rows that are assigned once.
          {
            j1 = static_cast<size_t>(rowSol[i]); //rowSol[i] is >= 0 here 
            min = -log(0);
            for (j = 0; j < dim; j++)  
              if (j != j1)
                if (assignCost(i, j - 1) - v[j] < min) 
                  min = assignCost(i, j - 1) - v[j - 1];
            v[j1] = v[j1] - min;
          }
        }
      }

      // Augmenting row reduction 
      short loopcnt = 0;           // do-loop to be done twice.
      do
      {
        loopcnt++;

        // scan all free rows.
        // in some cases, a free row may be replaced with another one to be scanned next.
        k = 0; 
        previousNumFree = numFree; 
        numFree = 0;             // start list of rows still free after augmenting row reduction.
        while (k < previousNumFree)
        {
          i = free[k]; 
          k++;

          // find minimum and second minimum reduced cost over columns.
          uMin = assignCost(i, 0) - v[0]; 
          j1 = 0; 
          uSubMin = static_cast<size_t>(-log(0));
          for (j = 1; j < dim; j++) 
          {
            h = assignCost(i, j) - v[j];
            if (h < uSubMin) {
              if (h >= uMin) 
              { 
                uSubMin = h; 
                j2 = j;
              }
              else 
              { 
                uSubMin = uMin; 
                uMin = h; 
                j2 = j1; 
                j1 = j;
              }
            }
          }

          i0 = colSol[j1];
          if (uMin < uSubMin) {
            // change the reduction of the minimum column to increase the minimum
            // reduced cost in the row to the subminimum.
            v[j1] = v[j1] - (uSubMin - uMin);
          } else {                  // minimum and subminimum equal.
            if (i0 >= 0)         // minimum column j1 is assigned.
            { 
              // swap columns j1 and j2, as j2 may be unassigned.
              j1 = j2; 
              i0 = colSol[j2];
            }
          }

          // (re-)assign i to j1, possibly de-assigning an i0.
          rowSol[i] = static_cast<int>(j1); 
          colSol[j1] = static_cast<int>(i);

          if (i0 >= 0) {          // minimum column j1 assigned earlier.
            if (uMin < uSubMin) {
              // put in current k, and go back to that k.
              // continue augmenting path i - j1 with i0.
              free[--k] = static_cast<size_t>(i0); 
            } else { 
              // no further augmenting reduction possible.
              // store i0 in list of free rows for next phase.
              free[numFree++] = static_cast<size_t>(i0); 
            }
          }
        }
      }
      while (loopcnt < 2);       // repeat once.

      // Augment solution for each free row.
      for (f = 0; f < numFree; f++) 
      {
        freeRow = free[f];       // start row of augmenting path.

        // Dijkstra shortest path algorithm.
        // runs until unassigned column added to shortest path tree.
        for (j = 0; j < dim; j++)  
        { 
          d[j] = assignCost(freeRow, j) - v[j]; 
          pred[j] = freeRow;
          colList[j] = j;        // init column list.
        }

        low = 0; // columns in 0..low-1 are ready, now none.
        up = 0;  // columns in low..up-1 are to be scanned for current minimum, now none.
        // columns in up..dim-1 are to be considered later to find new minimum, 
        // at this stage the list simply contains all columns 
        unassignedFound = false;
        do
        {
          if (up == low)         // no more columns to be scanned for current minimum.
          {
            last = low - 1; 

            // scan columns for up..dim-1 to find all indices for which new minimum occurs.
            // store these indices between low..up-1 (increasing up). 
            min = d[colList[up++]]; 
            for (k = up; k < dim; k++) 
            {
              j = colList[k]; 
              h = d[j];
              if (h <= min)
              {
                if (h < min)     // new minimum.
                { 
                  up = low;      // restart list at index low.
                  min = h;
                }
                // new index with same minimum, put on undex up, and extend list.
                colList[k] = colList[up]; 
                colList[up++] = j; 
              }
            }

            // check if any of the minimum columns happens to be unassigned.
            // if so, we have an augmenting path right away.
            for (k = low; k < up; k++) { 
              if (colSol[colList[k]] < 0) 
              {
                endOfPath = colList[k];
                unassignedFound = true;
                break;
              }
            }
          }

          if (!unassignedFound) 
          {
            // update 'distances' between freerow and all unscanned columns, via next scanned column.
            j1 = colList[low]; 
            low++; 
            i = static_cast<size_t>(colSol[j1]); 
            h = assignCost(i, j1) - v[j1] - min;

            for (k = up; k < dim; k++) 
            {
              j = colList[k]; 
              v2 = assignCost(i, j) - v[j] - h;
              if (v2 < d[j])
              {
                pred[j] = i;
                if (v2 == min) {  // new column found at same minimum value
                  if (colSol[j] < 0) 
                  {
                    // if unassigned, shortest augmenting path is complete.
                    endOfPath = j;
                    unassignedFound = true;
                    break;
                  }
                  // else add to list to be scanned right away.
                  else 
                  { 
                    colList[k] = colList[up]; 
                    colList[up++] = j; 
                  }
                }
                d[j] = v2;
              }
            }
          } 
        }
        while (!unassignedFound);

        // update column prices.
        for (k = 0; k <= last; k++)  
        { 
          j1 = colList[k]; 
          v[j1] = v[j1] + d[j1] - min;
        }

        // reset row and column assignments along the alternating path.
        do
        {
          i = pred[endOfPath]; 
          colSol[endOfPath] = static_cast<int>(i); 
          j1 = endOfPath; 
          endOfPath = static_cast<size_t>(rowSol[i]); 
          rowSol[i] = static_cast<int>(j1);
        }
        while (i != freeRow);
      }

      // calculate optimal cost.
      Scalar lapCost = 0;
      for (i = 0; i < dim; i++)  
      {
        j = static_cast<size_t>(rowSol[i]);
        u[i] = assignCost(i, j) - v[j];
        lapCost = lapCost + assignCost(i, j); 
      }

      return lapCost;
    }

  };

} // end of namespace bpp.

#endif  // _MATRIXTOOLS_H_
