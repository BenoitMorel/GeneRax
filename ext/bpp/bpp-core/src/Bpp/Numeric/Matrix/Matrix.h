//
// File: Matrix.h
// Authors: Julien Dutheil
//          Sylvain Gaillard
// Created on: Tue Apr 07 11:58 2004
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

   This software is a computer program whose purpose is to provide classes
   for numerical calculus.

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


#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <vector>
#include "../../Clonable.h"
#include "../NumConstants.h"
#include "../NumTools.h"
#include "../VectorExceptions.h"
#include <iostream>

namespace bpp
{
/**
 * @brief The matrix template interface.
 */
template<class Scalar>
class Matrix :
  public Clonable
{
public:
  Matrix() {}
  virtual ~Matrix() {}

public:

  /**
   * @return \f$m_{i,j}\f$.
   * @param i row index.
   * @param j column index.
   */
  virtual const Scalar& operator()(size_t i, size_t j) const = 0;
  
  /**
   * @return \f$m_{i,j}\f$.
   * @param i row index.
   * @param j column index.
   */
  virtual Scalar& operator()(size_t i, size_t j) = 0;

  virtual bool equals(const Matrix& m, double threshold = NumConstants::TINY())
  {
    if (m.getNumberOfRows() != getNumberOfRows() || m.getNumberOfColumns() != getNumberOfColumns())
      return false;
    for (size_t i = 0; i < getNumberOfRows(); i++)
    {
      for (size_t j = 0; j < getNumberOfColumns(); j++)
      {
        if (NumTools::abs<double>(static_cast<double>(operator()(i, j)) - static_cast<double>(m(i, j))) > threshold) return false;
      }
    }
    return true;
  }
  /**
   * @return The number of rows.
   */
  virtual size_t getNumberOfRows() const = 0;
  /**
   * @return The number of columns.
   */
  virtual size_t getNumberOfColumns() const = 0;
  /**
   * @return the row at position i as a vector.
   * @param i The index of the row.
   */
  virtual std::vector<Scalar> row(size_t i) const = 0;
  /**
   * @return the column at position j as a vector.
   * @param j The index of the column.
   */
  virtual std::vector<Scalar> col(size_t j) const = 0;
  /**
   * @brief Resize the matrix.
   *
   * @param nRows The new number of rows.
   * @param nCols The new number of columns.
   */
  virtual void resize(size_t nRows, size_t nCols) = 0;
};


/**
 * @brief Matrix storage by row.
 *
 * This matrix is a vector of vector of Scalar.
 * Row access is in \f$O(1)\f$ while column access is in \f$O(nRow)\f$.
 */
template<class Scalar>
class RowMatrix :
  public Matrix<Scalar>
{
protected:
  std::vector< std::vector<Scalar> > m_;

public:
  RowMatrix() : m_() {}

  RowMatrix(size_t nRow, size_t nCol) : m_(nRow)
  {
    for (size_t i = 0; i < nRow; i++)
    {
      m_[i].resize(nCol);
    }
  }

  RowMatrix(const Matrix<Scalar>& m) : m_(m.getNumberOfRows())
  {
    size_t nr = m.getNumberOfRows();
    size_t nc = m.getNumberOfColumns();
    for (size_t i = 0; i < nr; i++)
    {
      m_[i].resize(nc);
      for (size_t j = 0; j < nc; j++)
      {
        m_[i][j] = m(i, j);
      }
    }
  }

  RowMatrix& operator=(const Matrix<Scalar>& m)
  {
    size_t nr = m.getNumberOfRows();
    m_.resize(nr);
    size_t nc = m.getNumberOfColumns();
    for (size_t i = 0; i < nr; i++)
    {
      m_[i].resize(nc);
      for (size_t j = 0; j < nc; j++)
      {
        m_[i][j] = m(i, j);
      }
    }
    return *this;
  }

//  virtual ~RowMatrix() {}

public:
  RowMatrix* clone() const { return new RowMatrix(*this); }

  const Scalar& operator()(size_t i, size_t j) const { return m_[i][j]; }

  Scalar& operator()(size_t i, size_t j) { return m_[i][j]; }

  size_t getNumberOfRows() const { return m_.size(); }

  size_t getNumberOfColumns() const { return m_.size() == 0 ? 0 : m_[0].size(); }

  std::vector<Scalar> row(size_t i) const
  {
    std::vector<Scalar> r(getNumberOfColumns());
    for (size_t j = 0; j < getNumberOfColumns(); j++) { r[j] = operator()(i, j); }
    return r;
  }

  const std::vector<Scalar>& getRow(size_t i) const
  {
    return m_[i];
  }

  std::vector<Scalar>& getRow(size_t i) 
  {
    return m_[i];
  }

  std::vector<Scalar> col(size_t j) const
  {
    std::vector<Scalar> c(getNumberOfRows());
    for (size_t i = 0; i < getNumberOfRows(); i++) { c[i] = operator()(i, j); }
    return c;
  }

  void resize(size_t nRows, size_t nCols)
  {
    m_.resize(nRows);
    for (size_t i = 0; i < nRows; i++)
    {
      m_[i].resize(nCols);
    }
  }

  void addRow(const std::vector<Scalar>& newRow) noexcept(false)
  {
    if (getNumberOfColumns()!=0 && newRow.size() != getNumberOfColumns())
      throw DimensionException("RowMatrix::addRow: invalid row dimension", newRow.size(), getNumberOfColumns());
    m_.push_back(newRow);
  }
};

/**
 * @brief Matrix storage by column.
 *
 * This matrix is a vector of vector of Scalar.
 * Column access is in \f$O(1)\f$ while row access is in \f$O(nCol)\f$.
 */
  template<class Scalar>
  class ColMatrix :
    public Matrix<Scalar>
  {
  private:
    std::vector< std::vector<Scalar> > m_;

  public:
    ColMatrix() : m_() {}

    ColMatrix(size_t nRow, size_t nCol) : m_(nCol)
    {
      for (size_t i = 0; i < nCol; i++)
      {
        m_[i].resize(nRow);
      }
    }

    ColMatrix(const Matrix<Scalar>& m) : m_(m.getNumberOfColumns())
    {
      size_t nr = m.getNumberOfRows();
      size_t nc = m.getNumberOfColumns();
      for (size_t i = 0; i < nc; i++)
      {
        m_[i].resize(nr);
        for (size_t j = 0; j < nr; j++)
        {
          m_[i][j] = m(j, i);
        }
      }
    }

    ColMatrix& operator=(const Matrix<Scalar>& m)
    {
      size_t nc = m.getNumberOfColumns();
      m_.resize(nc);
      size_t nr = m.getNumberOfRows();
      for (size_t i = 0; i < nc; i++)
      {
        m_[i].resize(nr);
        for (size_t j = 0; j < nr; j++)
        {
          m_[i][j] = m(j, i);
        }
      }
      return *this;
    }

    virtual ~ColMatrix() {}

  public:
    ColMatrix* clone() const { return new ColMatrix(*this); }

    const Scalar& operator()(size_t i, size_t j) const { return m_[j][i]; }

    Scalar& operator()(size_t i, size_t j) { return m_[j][i]; }

    size_t getNumberOfColumns() const { return m_.size(); }

    size_t getNumberOfRows() const { return m_.size() == 0 ? 0 : m_[0].size(); }

    std::vector<Scalar> row(size_t i) const
    {
      std::vector<Scalar> r(getNumberOfColumns());
      for (size_t j = 0; j < getNumberOfColumns(); j++) { r[j] = operator()(i, j); }
      return r;
    }

    const std::vector<Scalar>& getCol(size_t i) const
    {
      return m_[i];
    }

    std::vector<Scalar> col(size_t j) const
    {
      std::vector<Scalar> c(getNumberOfRows());
      for (size_t i = 0; i < getNumberOfRows(); i++) { c[i] = operator()(i, j); }
      return c;
    }

    void resize(size_t nRows, size_t nCols)
    {
      m_.resize(nCols);
      for (size_t i = 0; i < nCols; i++)
      {
        m_[i].resize(nRows);
      }
    }

    void addCol(const std::vector<Scalar>& newCol) noexcept(false)
    {
      if (getNumberOfRows()!=0 && newCol.size() != getNumberOfRows())
        throw DimensionException("ColMatrix::addCol: invalid column dimension", newCol.size(), getNumberOfRows());
      m_.push_back(newCol);
    }
  };

/**
 * @brief Matrix storage in one vector.
 *
 * This Matrix is a simple vector of Scalar of size n x m.
 * Element access is in \f$O(1)\f$ but resizing the matrix while keeping the
 * old values is in \f$O(nm)\f$.
 *
 * Basic usage:
 * @code
 * LinearMatrix<int> m(3, 2); // Create a 3x2 matrix of int
 * m(1, 2) = 5; // Set the value of element at row = 1, col = 2 to 5
 * int x = m(0, 1); // Get the value of element at row = 0, col = 1;
 * @endcode
 *
 * @author Sylvain Gaillard
 */
template<class Scalar>
class LinearMatrix :
  public Matrix<Scalar>
{
private:
  std::vector<Scalar> m_;
  size_t rows_;
  size_t cols_;

public:
  /**
   * @brief Build a 0 x 0 matrix.
   */
  LinearMatrix() : m_(),
    rows_(0),
    cols_(0) { resize_(0, 0); }

  /**
   * @brief build a nRow x nCol matrix.
   */
  LinearMatrix(size_t nRow, size_t nCol) : m_(),
    rows_(nRow),
    cols_(nCol) { resize_(nRow, nCol); }

  LinearMatrix(const Matrix<Scalar>& m) : m_(m.getNumberOfRows() * m.getNumberOfColumns())
  {
    size_t nr = m.getNumberOfRows();
    size_t nc = m.getNumberOfColumns();
    for (size_t i = 0; i < nr; i++)
    {
      for (size_t j = 0; j < nc; j++)
      {
        m_[i * cols_ + j] = m(i, j);
      }
    }
  }

  LinearMatrix& operator=(const Matrix<Scalar>& m)
  {
    size_t nr = m.getNumberOfRows();
    size_t nc = m.getNumberOfColumns();
    m_.resize(nr * nc);
    for (size_t i = 0; i < nr; i++)
    {
      m_[i].resize(nc);
      for (size_t j = 0; j < nc; j++)
      {
        m_[i * cols_ + j] = m(i, j);
      }
    }
    return *this;
  }

  /**
   * @brief Destructor.
   */
  virtual ~LinearMatrix() {}

public:
  LinearMatrix* clone() const { return new LinearMatrix(*this); }

  const Scalar& operator()(size_t i, size_t j) const { return m_[i * cols_ + j]; }

  Scalar& operator()(size_t i, size_t j) { return m_[i * cols_ + j]; }

  size_t getNumberOfRows() const { return rows_; }

  size_t getNumberOfColumns() const { return cols_; }

  std::vector<Scalar> row(size_t i) const
  {
    std::vector<Scalar> r(getNumberOfColumns());
    for (size_t j = 0; j < getNumberOfColumns(); j++)
    {
      r[j] = operator()(i, j);
    }
    return r;
  }

  std::vector<Scalar> col(size_t j) const
  {
    std::vector<Scalar> c(getNumberOfRows());
    for (size_t i = 0; i < getNumberOfRows(); i++)
    {
      c[i] = operator()(i, j);
    }
    return c;
  }

  /**
   * @copydoc Matrix::resize
   *
   * This method resize the matrix keeping old data in place.
   * @see LinearMatrix::resize(size_t nRow, size_t nCol, bool keepValues)
   */
  void resize(size_t nRows, size_t nCols)
  {
    resize(nRows, nCols, true);
  }

  /**
   * @brief Resize the matrix.
   *
   * This task may be memory consumming if keepValues is true because it use
   * a copy of the input matrix to keep trace of the values.
   *
   * @param nRows the new number of rows
   * @param nCols the new number of columns
   * @param keepValues if old values must be kept in the resized matrix.
   * If keepValues = false, old values are still in the matrix but not at
   * the same positions. For instance:
   * @code
   * LinearMatrix<int> m(3, 2);
   * for (size_t i = 0 ; i < m.getNumberOfRows() ; i++) {
   *   for (size_t j = 0 ; j < m.getNumberOfColumns() ; j++) {
   *     m(i, j) = i * m.nCols() + j + 1;
   *   }
   * }
   * MatrixTools::print(m);
   * // 3x2
   * // [
   * // [1, 2]
   * // [3, 4]
   * // [5, 6]
   * // ]
   * LinearMatrix<int> m2 = m;
   * m2.resize(2, 4, false); // resize the matrix with keepValues = false
   * MatrixTools::print(m2);
   * // 2x4
   * // [
   * // [1, 2, 3, 4]
   * // [5, 6, 0, 0]
   * // ]
   * LinearMatrix<int> m3 = m;
   * m3.resize(2, 4, true); // resize the matrix with keepValues = true
   * MatrixTools::print(m3);
   * // 2x4
   * // [
   * // [1, 2, 0, 0]
   * // [3, 4, 0, 0]
   * // ]
   * @endcode
   */
  void resize(size_t nRows, size_t nCols, bool keepValues)
  {
    LinearMatrix<Scalar> tmpM;
    if (keepValues)
      tmpM = *this;
    resize_(nRows, nCols);
    if (keepValues)
    {
      for (size_t i = 0; i < nRows; i++)
      {
        for (size_t j = 0; j < nCols; j++)
        {
          if (i < tmpM.getNumberOfRows() && j < tmpM.getNumberOfColumns())
          {
            operator()(i, j) = tmpM(i, j);
          }
          else
          {
            operator()(i, j) = 0;
          }
        }
      }
    }
  }

private:
  /**
   * @brief Internal basic resize fonctionnalities.
   */
  void resize_(size_t nRows, size_t nCols)
  {
    m_.resize(nRows * nCols);
    rows_ = nRows;
    cols_ = nCols;
  }
};

template<class Scalar>
bool operator==(const Matrix<Scalar>& m1, const Matrix<Scalar>& m2)
{
  if (m1.getNumberOfRows() != m2.getNumberOfRows() || m1.getNumberOfColumns() != m2.getNumberOfColumns())
    return false;
  for (size_t i = 0; i < m1.getNumberOfRows(); i++)
  {
    for (size_t j = 0; j < m1.getNumberOfColumns(); j++)
    {
      if (m1(i, j) != m2(i, j))
        return false;
    }
  }
  return true;
}
} // end of namespace bpp.

#endif // _MATRIX_H_

