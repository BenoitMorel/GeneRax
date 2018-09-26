//
// File: VectorTools.h
// Created by: Julien Dutheil
// Created on: Fri Mar 14 14:16:32 2003
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

#ifndef _VECTORTOOLS_H_
#define _VECTORTOOLS_H_

#include "VectorExceptions.h"
#include "NumTools.h"
#include "AdaptiveKernelDensityEstimation.h"
#include "Matrix/Matrix.h"
#include "../Io/OutputStream.h"
#include "../App/ApplicationTools.h"

// From the STL:
#include <vector>
#include <map>
#include <cmath>
#include <algorithm>
#include <complex>

namespace bpp
{
typedef std::vector<std::complex<double> > Vcomplex;
typedef std::vector<Vcomplex> VVcomplex;
typedef std::vector<VVcomplex> VVVcomplex;

typedef std::vector<std::complex<long double> > Vlcomplex;
typedef std::vector<Vlcomplex> VVlcomplex;
typedef std::vector<VVlcomplex> VVVlcomplex;

typedef std::vector<double> Vdouble;
typedef std::vector<Vdouble> VVdouble;
typedef std::vector<VVdouble> VVVdouble;
typedef std::vector<VVVdouble> VVVVdouble;

typedef std::vector<long double> Vldouble;
typedef std::vector<Vldouble> VVldouble;
typedef std::vector<VVldouble> VVVldouble;
typedef std::vector<VVVldouble> VVVVldouble;

typedef std::vector<int> Vint;
typedef std::vector<Vint> VVint;
typedef std::vector<VVint> VVVint;
typedef std::vector<VVVint> VVVVint;

typedef std::vector<unsigned int> Vuint;
typedef std::vector<Vint> VVuint;
typedef std::vector<VVint> VVVuint;
typedef std::vector<VVVint> VVVVuint;

/**
 * @name Element-wise operations.
 * @{
 */

template<class T>
std::vector<T>  operator+(const std::vector<T>& v1, const std::vector<T>& v2)
noexcept(false)
{
  size_t size;
  if (v1.size() != v2.size())
  {
    throw DimensionException("VectorTools::operator+", v1.size(), v2.size());
  }
  else
  {
    size = v1.size();
  }
  std::vector<T> result(size);
  for (size_t i = 0; i < size; i++)
  {
    result[i] = v1[i] + v2[i];
  }
  return result;
}

template<class T>
std::vector<T> operator-(const std::vector<T>& v1, const std::vector<T>& v2)
noexcept(false)
{
  size_t size;
  if (v1.size() != v2.size())
  {
    throw DimensionException("VectorTools::operator-", v1.size(), v2.size());
  }
  else
  {
    size = v1.size();
  }
  std::vector<T> result(size);
  for (size_t i = 0; i < size; i++)
  {
    result[i] = v1[i] - v2[i];
  }
  return result;
}

template<class T>
std::vector<T> operator*(const std::vector<T>& v1, const std::vector<T>& v2)
noexcept(false)
{
  size_t size;
  if (v1.size() != v2.size())
  {
    throw DimensionException("VectorTools::operator*", v1.size(), v2.size());
  }
  else
  {
    size = v1.size();
  }
  std::vector<T> result(size);
  for (size_t i = 0; i < size; i++)
  {
    result[i] = v1[i] * v2[i];
  }
  return result;
}

template<class T>
std::vector<T> operator/(const std::vector<T>& v1, const std::vector<T>& v2)
noexcept(false)
{
  size_t size;
  if (v1.size() != v2.size())
  {
    throw DimensionException("VectorTools::operator/", v1.size(), v2.size());
  }
  else
  {
    size = v1.size();
  }
  std::vector<T> result(size);
  for (size_t i = 0; i < size; i++)
  {
    result[i] = v1[i] / v2[i];
  }
  return result;
}


template<class T, class C>
std::vector<T> operator+(const std::vector<T>& v1, const C& c)
{
  std::vector<T> result(v1.size());
  for (size_t i = 0; i < result.size(); i++)
  {
    result[i] = v1[i] + c;
  }
  return result;
}
template<class T, class C>
std::vector<T> operator+(const C& c, const std::vector<T>& v1)
{
  std::vector<T> result(v1.size());
  for (size_t i = 0; i < result.size(); i++)
  {
    result[i] = c + v1[i];
  }
  return result;
}

template<class T, class C>
std::vector<T> operator-(const std::vector<T>& v1, const C& c)
{
  std::vector<T> result(v1.size());
  for (size_t i = 0; i < result.size(); i++)
  {
    result[i] = v1[i] - c;
  }
  return result;
}
template<class T, class C>
std::vector<T> operator-(const C& c, const std::vector<T>& v1)
{
  std::vector<T> result(v1.size());
  for (size_t i = 0; i < result.size(); i++)
  {
    result[i] = c - v1[i];
  }
  return result;
}

template<class T, class C>
std::vector<T> operator*(const std::vector<T>& v1, const C& c)
{
  std::vector<T> result(v1.size());
  for (size_t i = 0; i < result.size(); i++)
  {
    result[i] = v1[i] * c;
  }
  return result;
}
template<class T, class C>
std::vector<T> operator*(const C& c, const std::vector<T>& v1)
{
  std::vector<T> result(v1.size());
  for (size_t i = 0; i < result.size(); i++)
  {
    result[i] = c * v1[i];
  }
  return result;
}

template<class T, class C>
std::vector<T> operator/(const std::vector<T>& v1, const C& c)
{
  std::vector<T> result(v1.size());
  for (size_t i = 0; i < result.size(); i++)
  {
    result[i] = v1[i] / c;
  }
  return result;
}
template<class T, class C>
std::vector<T> operator/(const C& c, const std::vector<T>& v1)
{
  std::vector<T> result(v1.size());
  for (size_t i = 0; i < result.size(); i++)
  {
    result[i] = c / v1[i];
  }
  return result;
}


template<class T>
void operator+=(std::vector<T>& v1, const std::vector<T>& v2)
noexcept(false)
{
  for (size_t i = 0; i < v1.size(); i++)
  {
    v1[i] += v2[i];
  }
}

template<class T>
void operator-=(std::vector<T>& v1, const std::vector<T>& v2)
noexcept(false)
{
  for (size_t i = 0; i < v1.size(); i++)
  {
    v1[i] -= v2[i];
  }
}

template<class T>
void operator*=(std::vector<T>& v1, const std::vector<T>& v2)
noexcept(false)
{
  for (size_t i = 0; i < v1.size(); i++)
  {
    v1[i] *= v2[i];
  }
}

template<class T>
void operator/=(std::vector<T>& v1, const std::vector<T>& v2)
noexcept(false)
{
  for (size_t i = 0; i < v1.size(); i++)
  {
    v1[i] /= v2[i];
  }
}

template<class T, class C>
void operator&=(std::vector<T>& v1, const C& c)
{
  for (size_t i = 0; i < v1.size(); i++)
  {
    v1[i] = c;
  }
}

template<class T, class C>
void operator+=(std::vector<T>& v1, const C& c)
{
  for (size_t i = 0; i < v1.size(); i++)
  {
    v1[i] += c;
  }
}

template<class T, class C>
void operator-=(std::vector<T>& v1, const C& c)
{
  for (size_t i = 0; i < v1.size(); i++)
  {
    v1[i] -= c;
  }
}

template<class T, class C>
void operator*=(std::vector<T>& v1, const C& c)
{
  for (size_t i = 0; i < v1.size(); i++)
  {
    v1[i] *= c;
  }
}

template<class T, class C>
void operator/=(std::vector<T>& v1, const C& c)
{
  for (size_t i = 0; i < v1.size(); i++)
  {
    v1[i] /= c;
  }
}
/** @} */

/******************************************************************************/

class VectorTools
{
public:
  VectorTools() {}
  virtual ~VectorTools() {}

public:
  /**
   * @name Matrix-like functions to resize arrays.
   *
   * @{
   */
  template<class T>
  static void resize2(VVdouble& vv, size_t n1, size_t n2)
  {
    vv.resize(n1);
    for (size_t i = 0; i < n1; i++) { vv[i].resize(n2); }
  }

  template<class T>
  static void resize3(VVVdouble& vvv, size_t n1, size_t n2, size_t n3)
  {
    vvv.resize(n1);
    for (size_t i = 0; i < n1; i++)
    {
      vvv[i].resize(n2);
      for (size_t j = 0; j < n2; j++)
      {
        vvv[i][j].resize(n3);
      }
    }
  }

  static void resize4(VVVVdouble& vvvv, size_t n1, size_t n2, size_t n3, size_t n4)
  {
    vvvv.resize(n1);
    for (size_t i = 0; i < n1; i++)
    {
      vvvv[i].resize(n2);
      for (size_t j = 0; j < n2; j++)
      {
        vvvv[i][j].resize(n3);
        for (size_t k = 0; k < n3; k++)
        {
          vvvv[i][j][k].resize(n4);
        }
      }
    }
  }
  /** @} */

  template<class T>
  static void fill(std::vector<T>& v, T value)
  {
    for (typename std::vector<T>::iterator it = v.begin(); it < v.end(); it++)
    {
      *it = value;
    }
  }

  /**
   * @brief Build a sequence std::vector.
   *
   * Build a std::vector from a value to another with a specified step.
   * This works for numerical values for which additions, subtractions and division
   * makes sens.
   *
   * @param from The begining.
   * @param to The end.
   * @param by The step.
   * @return A std::vector containing the sequence.
   */
  template<class T>
  static std::vector<T> seq(T from, T to, T by)
  {
    std::vector<T> v;
    if (from < to)
    {
      // for (T i = from ; i <= to ; i += by) {           // Not good for double, 'to'
      for (T i = from; i <= to + (by / 100); i += by)
      { // must be a little bit larger
        v.push_back(i);
      }
    }
    else
    {
      for (T i = from; i >= to - (by / 100); i -= by)
      {
        v.push_back(i);
      }
    }
    return v;
  }

  /**
   * @brief Send the position of the first occurence of 'which'.
   *
   * Comparisons are performed using the == operator.
   * Maximum complexity: O(v.size()).
   *
   * @param v The std::vector to search.
   * @param which The element to search.
   * @return The position of which in v.
   */
  template<class T>
  static size_t which(const std::vector<T>& v, const T& which)
  noexcept(false)
  {
    for (size_t i = 0; i < v.size(); i++)
    {
      if (v[i] == which) return i;
    }
    throw ElementNotFoundException<T>("VectorTools::which.", &v, &which);
  }

  /**
   * @brief Send the positions of all occurences of 'which'.
   *
   * Comparisons are performed using the == operator.
   * Complexity: O(v.size()).
   *
   * @param v The std::vector to search.
   * @param which The element to search.
   * @return A std::vector containing the positions of which in v.
   */
  template<class T>
  static std::vector<size_t> whichAll(const std::vector<T>& v, const T& which)
  noexcept(false)
  {
    std::vector<size_t> w;
    for (size_t i = 0; i < v.size(); i++)
    {
      if (v[i] == which) w.push_back(i);
    }
    if (w.size())
      return w;
    throw ElementNotFoundException<T>("VectorTools::whichAll.", &v, &which);
  }

  /**
   * @brief Send a new std::vector with unique elements.
   *
   * The input std::vector is copied, and the copy is sorted using QuickSort algorithm.
   * A one-pass loop then look for duplicates and copy unique element to a result std::vector.
   * The output std::vector is hence sorted.
   *
   * If v is empty, it is passed 'as is' in return (after being copied).
   *
   * @param v the std::vector to parse.
   */
  template<class T>
  static std::vector<T> unique(const std::vector<T>& v)
  {
    if (v.size() == 0) return v;
    std::vector<T> sortedV(v.begin(), v.end());
    sort(sortedV.begin(), sortedV.end());
    std::vector<T> uniq;
    uniq.push_back(sortedV[0]);
    for (size_t i = 1; i < sortedV.size(); i++)
    {
      if (sortedV[i] != sortedV[i - 1]) uniq.push_back(sortedV[i]);
    }
    return uniq;
  }

  /**
   * @brief Tell if the std::vector as unique elements.
   *
   * The input std::vector is copied, and the copy is sorted using QuickSort algorithm.
   * A one-pass loop then look for duplicates.
   *
   * If v is empty, the method returns 'true'.
   *
   * @param v the std::vector to parse.
   */
  template<class T>
  static bool isUnique(const std::vector<T>& v)
  {
    if (v.size() == 0) return true;
    std::vector<T> sortedV(v.begin(), v.end());
    sort(sortedV.begin(), sortedV.end());
    for (size_t i = 1; i < sortedV.size(); i++)
    {
      if (sortedV[i] == sortedV[i - 1]) return false;
    }
    return true;
  }

  /**
   * @author Laurent Gueguen
   * @param v1 the std::vector of elements,
   * @param v2 the std::vector of the selected positions
   * @return the std::vector of the selected elements, in the order of the
   *  required positions
   */
  template<class T>
  static std::vector<T> extract(const std::vector<T>& v1, const std::vector<int>& v2)
  {
    std::vector<T> v(v2.size());
    for (size_t i = 0; i < v2.size(); i++)
    {
      v[i] = v1[v2[i]];
    }
    return v;
  }

  /**
   * @brief Count each element of a std::vector.
   *
   * @return A map with keys = unique std::vector values and values = count for each std::vector value.
   * @param v the std::vector to parse.
   */
  template<class T>
  static std::map<T, size_t> countValues(const std::vector<T>& v)
  {
    std::map<T, size_t> c;
    for (size_t i = 0; i < v.size(); i++)
    {
      c[v[i]]++;
    }
    return c;
  }

  /**
   * @brief Get the break points for a given number of classes.
   *
   * Given a std::vector of values, return the values that cut the range of values
   * in a given number of classes.
   *
   * @param v The std::vector to parse.
   * @param n The expected number of classes.
   * @return a std::vector of size = n + 1 containing the breaking points.
   */
  static std::vector<double> breaks(const std::vector<double>& v, unsigned int n);

  /**
   * @brief Get the optimal class number following Scott's method.
   *
   * Use Scott's (1979) method to compute the optimal class number for histogram.
   *
   * Scott, D.W. (1979) On optimal and data-based histograms. Biometrika, 66, 605¿610.
   *
   * @param v The std::vector to parse.
   * @return The number of classes.
   */
  template<class T>
  static size_t nclassScott(const std::vector<T>& v)
  {
    std::vector<T> r1 = VectorTools::range(v);
    T r = r1[1] - r1[0];
    double n = v.size();
    double h = 3.5 * VectorTools::sd<T, double>(v) * std::pow(n, -1. / 3);
    return (size_t) ceil(r / h);
  }

  /**
   * @return The product of all elements in a std::vector.
   * @param v1 A std::vector.
   */
  template<class T>
  static T prod(const std::vector<T>& v1)
  {
    T p = 1;
    for (size_t i = 0; i < v1.size(); i++) { p *= v1[i]; }
    return p;
  }

  /**
   * @return The cumulative product of all elements in a std::vector.
   * @param v1 A std::vector.
   * @author Julien Dutheil
   */
  template<class T>
  static std::vector<T> cumProd(const std::vector<T>& v1)
  {
    std::vector<T> p(v1.size());
    if (v1.size() == 0) return p;
    p[0] = v1[0];
    for (size_t i = 1; i < v1.size(); i++) { p[i] = v1[i] * p[i - 1]; }
    return p;
  }

  /**
   * @return The sum of all elements in a std::vector.
   * @param v1 A std::vector.
   */
  template<class T>
  static T sum(const std::vector<T>& v1)
  {
    T s = 0;
    for (size_t i = 0; i < v1.size(); i++) { s += v1[i]; }
    return s;
  }

  /**
   * @return The cumulative sum of all elements in a std::vector.
   * @param v1 A std::vector.
   * @author Julien Dutheil
   */
  template<class T>
  static std::vector<T> cumSum(const std::vector<T>& v1)
  {
    std::vector<T> s(v1.size());
    if (v1.size() == 0) return s;
    s[0] = v1[0];
    for (size_t i = 1; i < v1.size(); i++) { s[i] = v1[i] + s[i - 1]; }
    return s;
  }

  /**
   * @author Laurent Gueguen
   * @return Log-normalize std::vector v1, ie add a constant to the elements of v
   *  such that @f$\sum_i(\exp(v_i)) = 1@f$.
   * @param v std::vector.
   */
  template<class T>
  static void logNorm(std::vector<T>& v)
  {
    T M = max(v);
    T x = std::exp(v[0] - M);
    for (size_t i = 1; i < v.size(); i++)
    {
      x += std::exp(v[i] - M);
    }
    v -= M + std::log(x);
  }

  /**
   * @author Laurent Gueguen
   * @return From std::vector v1, return @f$\log(\sum_i(\exp(v1_i)))@f$.
   * @param v1 a std::vector.
   */
  template<class T>
  static T logSumExp(const std::vector<T>& v1)
  {
    T M = max(v1);
    if (std::isinf(M))
      return M;
    
    T x = std::exp(v1[0] - M);
    for (size_t i = 1; i < v1.size(); i++)
    {
      x += std::exp(v1[i] - M);
    }
    return std::log(x) + M;
  }

  /**
   * @author Laurent Gueguen
   * @return From std::vector v1, return @f$\log(\sum_i(v2_i * \exp(v1_i)))@f$.
   * @param v1 a std::vector.
   * @param v2 another std::vector.
   */
  template<class T>
  static T logSumExp(const std::vector<T>& v1, const std::vector<T>& v2)
  {
    size_t size;
    if (v1.size() != v2.size())
      throw DimensionException("VectorTools::logsumexp", v1.size(), v2.size());
    else
      size = v1.size();

    T M = max(v1);
    if (std::isinf(M))
      throw BadNumberException("VectorTools::logSumExp", M);
    
    T x = v2[0] * std::exp(v1[0] - M);
    for (size_t i = 1; i < size; i++)
    {
      x += v2[i] * std::exp(v1[i] - M);
    }
    return std::log(x) + M;
  }

  /**
   * @author Laurent Gueguen
   * @return From std::vector v1, return @f$\log(\textrm{mean}_i(\exp(v1_i)))@f$.
   * @param v1 a std::vector.
   */
  template<class T>
  static T logMeanExp(const std::vector<T>& v1)
  {
    T M = max(v1);
    if (std::isinf(M))
      return M;

    T x = std::exp(v1[0] - M);
    for (size_t i = 1; i < v1.size(); i++)
    {
      x += std::exp(v1[i] - M);
    }
    return std::log(x) + M - std::log(v1.size());
  }


  /**
   * @author Laurent Gueguen
   * @return From std::vector v1, return @f$\sum_i(\exp(v1_i))@f$.
   * @param v1 a std::vector.
   */
  template<class T>
  static T sumExp(const std::vector<T>& v1)
  {
    T M = max(v1);
    if (std::isinf(M))
      return (M<0?0:M);

    T x = std::exp(v1[0] - M);
    for (size_t i = 1; i < v1.size(); i++)
    {
      x += std::exp(v1[i] - M);
    }
    return x * std::exp(M);
  }

  /**
   * @author Laurent Gueguen
   * @return From std::vector v1, return @f$\sum_i(v2_i * \exp(v1_i))@f$.
   * @param v1 a std::vector.
   * @param v2 another std::vector.
   */
  template<class T>
  static T sumExp(const std::vector<T>& v1, const std::vector<T>& v2)
  {
    size_t size;
    if (v1.size() != v2.size())
      throw DimensionException("VectorTools::sumExp", v1.size(), v2.size());
    else
      size = v1.size();

    T M = max(v1);
    if (std::isinf(M))
      throw BadNumberException("VectorTools::sumExp", M);

    T x = v2[0] * std::exp(v1[0] - M);
    for (size_t i = 1; i < size; i++)
    {
      x += v2[i] * std::exp(v1[i] - M);
    }
    return x * std::exp(M);
  }

  /**
   * @name These methods apply the corresponding function to each element
   * and return the result in a new std::vector.
   *
   * @{
   */
  
  template<class T>
  static std::vector<T> log(const std::vector<T>& v1)
  {
    std::vector<T> v2(v1.size());
    for (size_t i = 0; i < v2.size(); i++) { v2[i] = VectorTools::log(v1[i]); }
    return v2;
  }

  template<class T>
  static double log(const T& x)
  {
    return std::log(x);
  }
  
  template<class T>
  static std::vector<double> log(const std::vector<T>& v1, double base)
  {
    std::vector<double> v2(v1.size());
    for (size_t i = 0; i < v2.size(); i++) { v2[i] = std::log(v1[i]) / std::log(base); }
    return v2;
  }

  template<class T>
  static double exp(const T& x)
  {
    return std::exp(x);
  }

  template<class T>
  static std::vector<T> exp(const std::vector<T>& v1)
  {
    std::vector<T> v2(v1.size());
    for (size_t i = 0; i < v2.size(); i++) { v2[i] = VectorTools::exp(v1[i]); }
    return v2;
  }

  template<class T>
  static std::vector<double> cos(const std::vector<T>& v1)
  {
    std::vector<double> v2(v1.size());
    for (size_t i = 0; i < v2.size(); i++) { v2[i] = std::cos(v1[i]); }
    return v2;
  }

  template<class T>
  static std::vector<double> sin(const std::vector<T>& v1)
  {
    std::vector<double> v2(v1.size());
    for (size_t i = 0; i < v2.size(); i++) { v2[i] = std::sin(v1[i]); }
    return v2;
  }

  template<class T>
  static std::vector<double> log10(const std::vector<T>& v1)
  {
    std::vector<double> v2(v1.size());
    for (size_t i = 0; i < v1.size(); i++) { v2[i] = std::log10(v1[i]); }
    return v2;
  }

  template<class T>
  static std::vector<T> fact(const std::vector<T>& v1)
  {
    std::vector<T> v2(v1.size());
    for (size_t i = 0; i < v1.size(); i++) { v2[i] = NumTools::fact<T>(v1[i]); }
    return v2;
  }

  template<class T>
  static std::vector<T> sqr(const std::vector<T>& v1)
  {
    std::vector<T> v2(v1.size());
    for (size_t i = 0; i < v1.size(); i++) { v2[i] = NumTools::sqr<T>(v1[i]); }
    return v2;
  }

  template<class T>
  static std::vector<T> pow(const std::vector<T>& v1, T& b)
  {
    std::vector<T> v2(v1.size());
    for (size_t i = 0; i < v1.size(); i++) { v2[i] = std::pow(v1[i], b); }
    return v2;
  }
  /** @} */

  /**
   * @brief Concatenate a std::vector after converting to string.
   *
   * @param v The std::vector to concatenate.
   * @param delim A string which is used to separate the values (default is " ").
   */
  template<class T>
  static std::string paste(const std::vector<T>& v, const std::string& delim = " ")
  {
    std::ostringstream out;
    out.precision(12);
    for (size_t i = 0; i < v.size(); i++)
    {
      out << v[i];
      if (i < v.size() - 1)
        out << delim;
    }
    return out.str();
  }

  /**
   * @brief Print a std::vector to a stream.
   * @param v1 A std::vector.
   * @param out A stream.
   * @param delim A string which is used to separate the values (default is " ").
   */
  template<class T>
  static void print(const std::vector<T>& v1, OutputStream& out = * ApplicationTools::message, const std::string& delim = " ")
  {
    for (size_t i = 0; i < v1.size(); i++)
    {
      out << v1[i];
      if (i < v1.size() - 1)
        out << delim;
    }
    out.endLine();
  }

  template<class T>
  static void printRange(const std::vector<T>& v1, OutputStream& out = * ApplicationTools::message, const std::string& delim = ",", const std::string& rangeDelim = "-")
  {
    size_t vs=v1.size();
    
    for (size_t i = 0; i < vs; i++)
    {
      out << v1[i];
      size_t j=i+1;
      
      while (j<vs)
        if (v1[j]==v1[j-1]+1)
          j++;
        else
          break;
      
      if (j>i+2)
        out << rangeDelim << v1[j-1];
      i=j-1;
      if (i<vs-1)
        out << delim;
    }
  }


  /**
   * @brief Print a std::vector to a stream in R format
   * @param v1 A std::vector.
   * @param variableName the variable name (default "x").
   * @param out A stream (default cout).
   */
  template<class T>
  static void printForR(const std::vector<T>& v1, std::string variableName = "x", std::ostream& out = std::cout)
  {
    out.precision(12);
    out << variableName << "<-c(";
    for (size_t i = 0; i < v1.size(); i++)
    {
      out << v1[i];
      if (i < v1.size() - 1)
        out << ", ";
    }
    out << ")" << std::endl;
  }

  /**
   * @return The scalar product of two std::vectors.
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   * @throw DimensionException If the two std::vector do not have the same length.
   */
  template<class InputType, class OutputType>
  static OutputType scalar(const std::vector<InputType>& v1, const std::vector<InputType>& v2)
  noexcept(false)
  {
    if (v1.size() != v2.size())
    {
      throw DimensionException("VectorTools::scalar", v1.size(), v2.size());
    }
    OutputType result = 0;
    for (size_t i = 0; i < v1.size(); i++)
    {
      result += v1[i] * v2[i];
    }
    return result;
  }
  /**
   * This dt product correspond to the dot product <v1,v2> in the space defined by
   * @f[
   * M =
   * \begin{pmatrix}
   * w_1 & \ldots & \\
   * \vdots & w_2  & \ldots\\
   *        & \vdots & \ddots\\
   * \end{pmatrix}
   * @f]
   * @return The "weighted" scalar product of two std::vectors.
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   * @param w A std::vector of weights.
   * @throw DimensionException If the two std::vector do not have the same length or do not match the length of the weights.
   */
  template<class InputType, class OutputType>
  static OutputType scalar(const std::vector<InputType>& v1, const std::vector<InputType>& v2, const std::vector<InputType>& w)
  noexcept(false)
  {
    if (v1.size() != w.size())
    {
      throw DimensionException("VectorTools::scalar", v1.size(), w.size());
    }
    if (v2.size() != w.size())
    {
      throw DimensionException("VectorTools::scalar", v2.size(), w.size());
    }
    OutputType result = 0;
    for (size_t i = 0; i < v1.size(); i++)
    {
      result += v1[i] * v2[i] * w[i];
    }
    return result;
  }

  /**
   * @return The scalar Kronecker product of two std::vectors.
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   * @throw DimensionException If the two std::vector do not have the same length.
   */
  template<class T>
  static std::vector<T> kroneckerMult(const std::vector<T>& v1, const std::vector<T>& v2)
  noexcept(false)
  {
    size_t n1 = v1.size();
    size_t n2 = v2.size();
    std::vector<T> v3(n1 * n2);
    for (size_t i = 0; i < n1; i++)
    {
      T v1i = v1[i];
      for (size_t j = 0; j < n2; j++)
      {
        v3[i * n2 + j] = v1i * v2[j];
      }
    }
    return v3;
  }

  /**
   * @return The norm of a std::vector (@f$\sqrt{\sum_i^n x_i^2}@f$).
   * @param v1 A std::vector.
   */
  template<class InputType, class OutputType>
  static OutputType norm(const std::vector<InputType>& v1)
  {
    OutputType result = 0;
    for (size_t i = 0; i < v1.size(); i++)
    {
      result += v1[i] * v1[i];
    }
    return sqrt(result);
  }

  /**
   * @return The "weighted" norm of a std::vector (@f$\sqrt{\sum_i^n x_i^2}@f$).
   * @param v1 A std::vector.
   * @param w A std::vector of weights.
   * @throw DimensionException If v1 and w do not have the same length.
   * @see scalar.
   */
  template<class InputType, class OutputType>
  static OutputType norm(const std::vector<InputType>& v1, const std::vector<InputType>& w)
  noexcept(false)
  {
    if (v1.size() != w.size())
    {
      throw DimensionException("VectorTools::norm", v1.size(), w.size());
    }
    OutputType result = 0;
    for (size_t i = 0; i < v1.size(); i++)
    {
      result += v1[i] * v1[i] * w[i];
    }
    return sqrt(result);
  }

  /**
   * @return The cosinus of the angle of two std::vectors.
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   * @throw DimensionException If the two std::vector do not have the same length.
   */
  template<class InputType, class OutputType>
  static OutputType cos(const std::vector<InputType>& v1, const std::vector<InputType>& v2)
  noexcept(false)
  {
    return scalar<InputType, OutputType>(v1, v2)
           / (norm<InputType, OutputType>(v1) * norm<InputType, OutputType>(v2));
  }

  /**
   * @return The weighted cosinus of the angle of two std::vectors.
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   * @param w A std::vector of weights.
   * @throw DimensionException If the two std::vector do not have the same length.
   */
  template<class InputType, class OutputType>
  static OutputType cos(const std::vector<InputType>& v1, const std::vector<InputType>& v2, const std::vector<InputType>& w)
  noexcept(false)
  {
    return scalar<InputType, OutputType>(v1, v2, w)
           / (norm<InputType, OutputType>(v1, w) * norm<InputType, OutputType>(v2, w));
  }

  /**
   * @name Extrema.
   *
   * @{
   */

  /**
   * @brief Template function to get the minimum value of a std::vector.
   *
   * The < operator must be defined for the specified class.
   *
   * @param v The input std::vector.
   * @return The minimum value in the std::vector.
   * @throw EmptyVectorException If the input std::vector is empty.
   */
  template<class T>
  static T min(const std::vector<T>& v) noexcept(false)
  {
    if (v.size() == 0) throw EmptyVectorException<T>("VectorTools::min()", &v);
    T mini = v[0];
    for (size_t i = 1; i < v.size(); i++)
    {
      if (v[i] < mini) mini = v[i];
    }
    return mini;
  }

  /**
   * @brief Template function to get the maximum value of a std::vector.
   *
   * The > operator must be defined for the specified class.
   *
   * @param v The input std::vector.
   * @return The maximum value in the std::vector.
   * @throw EmptyVectorException If the input std::vector is empty.
   */
  template<class T>
  static T max(const std::vector<T>& v) noexcept(false)
  {
    if (v.size() == 0) throw EmptyVectorException<T>("VectorTools::max()", &v);
    T maxi = v[0];
    for (size_t i = 1; i < v.size(); i++)
    {
      if (v[i] > maxi) maxi = v[i];
    }
    return maxi;
  }

  /**
   * @brief Template function to get the index of the maximum value of a std::vector.
   *
   * The > operator must be defined for the specified class.
   * The position sent is the first one matching the maximum value.
   *
   * @param v The input std::vector.
   * @return The position of the maximum value in the std::vector.
   * @throw EmptyVectorException If the input std::vector is empty.
   */
  template<class T>
  static size_t whichMax(const std::vector<T>& v) noexcept(false)
  {
    if (v.size() == 0) throw EmptyVectorException<T>("VectorFuntions::whichMax()", &v);
    T maxi = v[0];
    size_t pos = 0;
    for (size_t i = 1; i < v.size(); i++)
    {
      if (v[i] > maxi)
      {
        maxi = v[i];
        pos = i;
      }
    }
    return pos;
  }

  /**
   * @brief Template function to get the index of the minimum value of a std::vector.
   *
   * The < operator must be defined for the specified class.
   * The position sent is the first one matching the maximum value.
   *
   * @param v The input std::vector.
   * @return The position of the minimum value in the std::vector.
   * @throw EmptyVectorException If the input std::vector is empty.
   */
  template<class T>
  static size_t whichMin(const std::vector<T>& v) noexcept(false)
  {
    if (v.size() == 0) throw EmptyVectorException<T>("VectorTools::whichMin()", &v);
    T mini = v[0];
    size_t pos = 0;
    for (size_t i = 1; i < v.size(); i++)
    {
      if (v[i] < mini)
      {
        mini = v[i];
        pos = i;
      }
    }
    return pos;
  }

  /**
   * @brief Template function to get the indices of the maximum value of a std::vector.
   *
   * The > and == operator must be defined for the specified class.
   * All positions matching the maximum value are returned.
   *
   * @param v The input std::vector.
   * @return A vector of positions matching the maximum value in the std::vector.
   * @throw EmptyVectorException If the input std::vector is empty.
   */
  template<class T>
  static std::vector<size_t> whichMaxAll(const std::vector<T>& v)
  noexcept(false)
  {
    if (v.size() == 0) throw EmptyVectorException<T>("VectorFuntions::whichMaxAll()", &v);
    T maxi = max(v);
    std::vector<size_t> pos;
    for (size_t i = 0; i < v.size(); i++)
    {
      if (v[i] == maxi)
      {
        pos.push_back(i);
      }
    }
    return pos;
  }

  /**
   * @brief Template function to get the indices of the minimum value of a std::vector.
   *
   * The < and == operator must be defined for the specified class.
   * All positions matching the maximum value are returned.
   *
   * @param v The input std::vector.
   * @return A vector of positions matching the minimum value in the std::vector.
   * @throw EmptyVectorException If the input std::vector is empty.
   */
  template<class T>
  static std::vector<size_t> whichMinAll(const std::vector<T>& v)
  noexcept(false)
  {
    if (v.size() == 0) throw EmptyVectorException<T>("VectorFuntions::whichMinAll()", &v);
    T mini = min(v);
    std::vector<size_t> pos;
    for (size_t i = 0; i < v.size(); i++)
    {
      if (v[i] == mini)
      {
        pos.push_back(i);
      }
    }
    return pos;
  }


  /**
   * @brief Template function to get both extrema of a std::vector.
   *
   * Both < and > operators must be defined for the specified class.
   *
   * @param v The input std::vector.
   * @return A std::vector of size 2 which values are min(v) and max(v).
   * throw EmptyVectorException If the input std::vector is empty.
   */
  template<class T>
  static std::vector<T> range(const std::vector<T>& v) noexcept(false)
  {
    if (v.size() == 0)
      throw EmptyVectorException<T>("VectorTools::range()", &v);
    std::vector<T> r(2);
    r[0] = r[1] = v[0];
    for (size_t i = 1; i < v.size(); i++)
    {
      if (v[i] < r[0]) r[0] = v[i];
      if (v[i] > r[1]) r[1] = v[i];
    }
    return r;
  }

private:
  template<class T> class order_Cmp_
  {
    const std::vector<T>& values_;

public:
    order_Cmp_(const std::vector<T>& v) : values_(v) {}
    bool operator()(size_t a, size_t b) { return values_[a] < values_[b]; }
  };

public:
  /**
   * @brief Template function to get the order of elements in the input vector.
   *
   * This function is equivalent to the R 'order' function. It will sort elements by increasing size, so that [0] in the resulting outptu vector is the indice of the minimum element in the input vector, [1] is the second minimum, etc.
   *
   * @param v The input std::vector.
   * @return Vector with sorted positions, according to input vector.
   * throw EmptyVectorException If the input std::vector is empty.
   */
  template<class T>
  static std::vector<size_t> order(const std::vector<T>& v) noexcept(false)
  {
    if (v.size() == 0)
      throw EmptyVectorException<T>("VectorTools::sort()", &v);
    // Private inner class:
    std::vector<size_t> index(v.size());
    for (size_t i = 0; i < index.size(); ++i)
    {
      index[i] = i;
    }
    sort(index.begin(), index.end(), order_Cmp_<T>(v));
    return index;
  }

  /** @} */

  /**
   * @brief Template function to get the absolute value of all elements of a std::vector.
   *
   *
   * @param v The input std::vector.
   * @return A vector with all absolute values.
   */
  template<class T>
  static std::vector<T> abs(const std::vector<T>& v) noexcept(false)
  {
    std::vector<T> vabs(v.size());
    for (size_t i = 1; i < v.size(); i++)
    {
      vabs[i] = std::abs(v[i]);
    }
    return vabs;
  }

  /**
   * @return The mean value of the std::vector.
   * @param v1 A std::vector.
   */
  template<class InputType, class OutputType>
  static OutputType mean(const std::vector<InputType>& v1)
  {
    return (OutputType)sum<InputType>(v1) / (OutputType)v1.size();
  }
  /**
   * @return The weighted mean value of the std::vector.
   * @param v1 A std::vector.
   * @param w A std::vector of weights.
   * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
   */
  template<class InputType, class OutputType>
  static OutputType mean(const std::vector<InputType>& v1, const std::vector<InputType>& w, bool normalizeWeights = true)
  {
    if (normalizeWeights)
    {
      std::vector<InputType> wn = w / sum(w);
      return scalar<InputType, OutputType>(v1, wn);
    }
    else
    {
      return scalar<InputType, OutputType>(v1, w);
    }
  }

  /**
   * @return The median value of the std::vector.
   * @param v1 A std::vector.
   */
  template<class InputType>
  static InputType median(std::vector<InputType>& v1)
  {
    InputType med = 0;
    if (v1.size() == 0) return med;
    if (v1.size() == 1) return v1[0];
    sort(v1.begin(), v1.end());
    size_t i = v1.size() / 2;
    if (v1.size() % 2 == 0)
    {
      // Vector size is pair
      med = double((v1[i - 1] + v1[i]) / 2);
    }
    else
    {
      // Vector size is impair
      med = v1[i];
    }
    return med;
  }

  /**
   * @brief Set the mean of a std::vector to be 0.
   *
   * @return A std::vector with mean 0.
   * @param v1 A std::vector.
   */
  template<class InputType, class OutputType>
  static std::vector<OutputType> center(const std::vector<InputType>& v1)
  {
    OutputType m = mean<InputType, OutputType>(v1);
    std::vector<OutputType> v(v1.size());
    for (size_t i = 0; i < v1.size(); i++)
    {
      v[i] = (OutputType)v1[i] - m;
    }
    return v;
  }
  /**
   * @brief Set the weighted mean of a std::vector to be 0.
   *
   * @return A std::vector with mean 0.
   * @param v1 A std::vector.
   * @param w A std::vector of weights.
   * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
   */
  template<class InputType, class OutputType>
  static std::vector<OutputType> center(const std::vector<InputType>& v1, const std::vector<InputType>& w, bool normalizeWeights = true)
  {
    OutputType m = mean<InputType, OutputType>(v1, w, normalizeWeights);
    std::vector<OutputType> v(v1.size());
    for (size_t i = 0; i < v1.size(); i++)
    {
      v[i] = (OutputType)v1[i] - m;
    }
    return v;
  }

  /**
   * @return The covariance of two std::vectors.
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   * @param unbiased Tell if an unbiased estimate must be computed.
   * @throw DimensionException If the two std::vector do not have the same length.
   */
  template<class InputType, class OutputType>
  static OutputType cov(const std::vector<InputType>& v1, const std::vector<InputType>& v2, bool unbiased = true)
  noexcept(false)
  {
    OutputType n = (OutputType)v1.size();
    OutputType x =  scalar<InputType, OutputType>(
      center<InputType, OutputType>(v1),
      center<InputType, OutputType>(v2)
      ) / n;
    if (unbiased) x = x * n / (n - 1);
    return x;
  }

  /**
   * @return The weighted covariance of two std::vectors.
   * To have a population estimate you have to multiply by \f$\frac{n}{n-1}\f$.
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   * @param w A std::vector of weights.
   * @param unbiased Tell if an unbiased estimate must be computed.
   * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
   * @throw DimensionException If the two std::vector do not have the same length.
   */
  template<class InputType, class OutputType>
  static OutputType cov(const std::vector<InputType>& v1, const std::vector<InputType>& v2, const std::vector<InputType>& w, bool unbiased = true, bool normalizeWeights = true)
  noexcept(false)
  {
    if (normalizeWeights)
    {
      std::vector<InputType> wn = w / sum(w);
      OutputType x = scalar<InputType, OutputType>(
        center<InputType, OutputType>(v1, wn, false),
        center<InputType, OutputType>(v2, wn, false),
        wn
        );
      if (unbiased)
      {
        x = x / (1 - sum(sqr<double>(wn)));
      }
      return x;
    }
    else
    {
      OutputType x = scalar<InputType, OutputType>(
        center<InputType, OutputType>(v1, w, false),
        center<InputType, OutputType>(v2, w, false),
        w
        );
      if (unbiased)
      {
        x = x / (1 - sum(sqr(w)));
      }
      return x;
    }
  }
  /**
   * @return The variance of the std::vector.
   * @param v1 The sample std::vector.
   * @param unbiased Tell if an unbiased estimate must be computed.
   */
  template<class InputType, class OutputType>
  static OutputType var(const std::vector<InputType>& v1, bool unbiased = true)
  {
    return cov<InputType, OutputType>(v1, v1, unbiased);
  }
  /**
   * @return The weighted variance of the std::vector.
   * @param v1 The sample std::vector.
   * @param w A std::vector of weights.
   * @param unbiased Tell if an unbiased estimate must be computed.
   * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
   * @throw DimensionException If v1 and w do not have the same length.
   */
  template<class InputType, class OutputType>
  static OutputType var(const std::vector<InputType>& v1, const std::vector<InputType>& w, bool unbiased = true, bool normalizeWeights = true)
  noexcept(false)
  {
    return cov<InputType, OutputType>(v1, v1, w, unbiased, normalizeWeights);
  }

  /**
   * @return The standard deviation of the std::vector.
   * @param v1 The sample std::vector.
   * @param unbiased Tell if an unbiased estimate must be computed.
   */
  template<class InputType, class OutputType>
  static OutputType sd(const std::vector<InputType>& v1, bool unbiased = true)
  {
    return sqrt(var<InputType, OutputType>(v1, unbiased));
  }

  /**
   * @return The weighted standard deviation of the std::vector.
   * @param v1 The sample std::vector.
   * @param w A std::vector of weights.
   * @param unbiased Tell if an unbiased estimate must be computed.
   * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
   * @throw DimensionException If v1 and w do not have the same length.
   */
  template<class InputType, class OutputType>
  static OutputType sd(const std::vector<InputType>& v1, const std::vector<InputType>& w, bool unbiased = true, bool normalizeWeights = true)
  noexcept(false)
  {
    return sqrt(var<InputType, OutputType>(v1, w, unbiased, normalizeWeights));
  }

  /**
   * @return The Pearson correlation coefficient of two std::vectors.
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   * @throw DimensionException If the two std::vector do not have the same length.
   */
  template<class InputType, class OutputType>
  static OutputType cor(const std::vector<InputType>& v1, const std::vector<InputType>& v2)
  noexcept(false)
  {
    return cov<InputType, OutputType>(v1, v2)
           / ( sd<InputType, OutputType>(v1) * sd<InputType, OutputType>(v2) );
  }

  /**
   * @return The weighted Pearson correlation coefficient of two std::vectors.
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   * @param w A std::vector of weights.
   * @param normalizeWeights Tell if weights should be normalized so that they sum to 1.
   * @throw DimensionException If the two std::vector do not have the same length.
   */
  template<class InputType, class OutputType>
  static OutputType cor(const std::vector<InputType>& v1, const std::vector<InputType>& v2, const std::vector<InputType>& w, bool normalizeWeights = true)
  noexcept(false)
  {
    if (normalizeWeights)
    {
      std::vector<InputType> wn = w / sum(w);
      return cov<InputType, OutputType>(v1, v2, wn, false, false)
             / ( sd<InputType, OutputType>(v1, wn, false, false) * sd<InputType, OutputType>(v2, wn, false, false) );
    }
    else
    {
      return cov<InputType, OutputType>(v1, v2, w, false, false)
             / ( sd<InputType, OutputType>(v1, w, false, false) * sd<InputType, OutputType>(v2, w, false, false) );
    }
  }

  /**
   * @return Compute the Shannon entropy indice of a vector of frequencies.
   *
   * The elements of the vector should contain frequencies
   * of each modality of the variable and sum to one. This verification is left to the user
   * though.
   *
   * @author Julien Dutheil
   * @see shannonDiscrete for computing the shannon indice of a sample, discrete version.
   * @see shannonContinuous for computing the shannon indice of a sample, continuous version.
   *
   * @param v The input std::vector.
   * @param base The base of the logarithm to use.
   */
  template<class InputType, class OutputType>
  static OutputType shannon(const std::vector<InputType>& v, double base = 2.7182818)
  {
    OutputType s = 0;
    for (size_t i = 0; i < v.size(); i++)
    {
      if (v[i] > 0) s += static_cast<OutputType>(v[i] * std::log(v[i]) / std::log(base));
    }
    return -s;
  }

  /**
   * @return Compute the Shannon entropy indice of a vector.
   *
   * This is the discrete version. The vector is supposed to be a finite sample from
   * a discrete distribution. The counts of each observed state are computed and used
   * to compute the discrete entropy.
   *
   * @author Julien Dutheil
   * @see shannonContinuous for the continuous version.
   * @see shannon for direct computation from frequencies.
   *
   * @param v The input std::vector.
   * @param base The base of the logarithm to use.
   */
  template<class InputType, class OutputType>
  static OutputType shannonDiscrete(const std::vector<InputType>& v, double base = 2.7182818)
  {
    std::map<InputType, double> counts;
    for (size_t i = 0; i < v.size(); i++)
    {
      counts[v[i]]++;
    }
    OutputType s = 0;
    double n = static_cast<double>(v.size());
    for (typename std::map<InputType, double>::iterator it = counts.begin(); it != counts.end(); it++)
    {
      s += static_cast<OutputType>((it->second / n) * std::log(it->second / n) / std::log(base));
    }
    return -s;
  }

  /**
   * @return Compute mutual information index from two samples from discrete variables.
   *
   * This is the discrete version. Each vector is supposed to be a finite sample from
   * a discrete distribution. The counts of each (joint) observed state are computed
   * and used to compute the discrete mutual information.
   *
   * @author Julien Dutheil
   * @see miContinuous.
   *
   * @param v1 The first input vector.
   * @param v2 The second input vector.
   * @param base The base of the logarithm to use.
   * @throw DimensionException if the two vectors do not have the same lengths.
   */
  template<class InputType, class OutputType>
  static OutputType miDiscrete(const std::vector<InputType>& v1, const std::vector<InputType>& v2, double base = 2.7182818)
  noexcept(false)
  {
    if (v1.size() != v2.size())
      throw DimensionException("VectorTools::miDiscrete. The two samples must have the same length.", v2.size(), v1.size());
    std::map<InputType, double> counts1;
    std::map<InputType, double> counts2;
    std::map<InputType, std::map<InputType, double> > counts12;
    for (size_t i = 0; i < v1.size(); i++)
    {
      counts1[v1[i]]++;
      counts2[v2[i]]++;
      counts12[v1[i]][v2[i]]++;
    }
    OutputType s = 0;
    double n = static_cast<double>(v1.size());
    for (typename std::map<InputType, std::map<InputType, double> >::iterator it1 = counts12.begin(); it1 != counts12.end(); it1++)
    {
      for (typename std::map<InputType, double>::iterator it2 = it1->second.begin(); it2 != it1->second.end(); it2++)
      {
        s += static_cast<OutputType>((it2->second / n) * std::log(it2->second * n / (counts1[it1->first] * counts2[it2->first])) / std::log(base));
      }
    }
    return s;
  }

  /**
   * @return Compute the Shannon entropy indice of a sample from a continuous variable.
   *
   * This is the continuous version. The vector is supposed to be a finite sample from
   * a continuous distribution. The density is of the distribution is estimated using
   * a kernel method, and is used to compute the continuous entropy.
   *
   * Reference: Ivan Kojadinovic (2004) _Computational Statistics & Data Analysis_, 46:269-294
   *
   * @author Julien Dutheil
   * @see shannon For the discrete version.
   *
   * @param v The input std::vector.
   * @param base The base of the logarithm to use.
   */
  template<class InputType, class OutputType>
  static OutputType shannonContinuous(const std::vector<InputType>& v, double base = 2.7182818)
  {
    LinearMatrix<InputType> m(1, v.size());
    for (size_t i = 0; i < v.size(); i++)
    {
      m(0, i) = v[i];
    }
    AdaptiveKernelDensityEstimation kd(m);
    OutputType s = 0;
    std::vector<double> x(1);
    for (size_t i = 0; i < v.size(); i++)
    {
      x[0] = static_cast<double>(v[i]);
      s += static_cast<OutputType>(std::log(kd.kDensity(x)) / std::log(base));
    }
    return -s / static_cast<double>(v.size());
  }

  /**
   * @return Compute mutual information index from two samples from continuous variables.
   *
   * This is the continuous version. Each vector is supposed to be a finite sample from
   * a continuous distribution. The density is of the distribution is estimated using
   * a kernel method, as well as the joint density, and are used to compute the continuous
   * mutual information.
   *
   * Reference: Ivan Kojadinovic (2004) _Computational Statistics & Data Analysis_, 46:269-294
   *
   * @author Julien Dutheil
   * @see miDiscrete For the discrete version.
   * @warning This function is experimental. Notably, bad kernel density estimates may lead to negative MI :(
   *
   * @param v1 The first input vector.
   * @param v2 The second input vector.
   * @param base The base of the logarithm to use.
   * @throw DimensionException if the two vectors do not have the same lengths.
   */
  template<class InputType, class OutputType>
  static OutputType miContinuous(const std::vector<InputType>& v1, const std::vector<InputType>& v2, double base = 2.7182818)
  noexcept(false)
  {
    if (v1.size() != v2.size())
      throw DimensionException("VectorTools::miContinuous. The two samples must have the same length.", v2.size(), v1.size());
    LinearMatrix<InputType> m1(1, v1.size());
    LinearMatrix<InputType> m2(1, v2.size());
    LinearMatrix<InputType> m12(2, v1.size());
    for (size_t i = 0; i < v1.size(); i++)
    {
      m1(0, i) = m12(0, i) = v1[i];
      m2(0, i) = m12(1, i) = v2[i];
    }
    AdaptiveKernelDensityEstimation kd1(m1);
    AdaptiveKernelDensityEstimation kd2(m2);
    AdaptiveKernelDensityEstimation kd12(m12);
    OutputType s = 0;
    std::vector<double> x1(1);
    std::vector<double> x2(1);
    std::vector<double> x12(2);
    for (size_t i = 0; i < v1.size(); i++)
    {
      x1[0] = x12[0] = static_cast<double>(v1[i]);
      x2[0] = x12[1] = static_cast<double>(v2[i]);
      s += static_cast<OutputType>(std::log(kd12.kDensity(x12) / (kd1.kDensity(x1) * kd2.kDensity(x2))) / std::log(base));
    }
    return s / static_cast<double>(v1.size());
  }

  /**
   * @return 'true' if the two std::vectors contains the same elements, whatever their order in the container.
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   */
  template<class T>
  static bool haveSameElements(const std::vector<T>& v1, const std::vector<T>& v2)
  {
    std::vector<T> u1(v1);
    std::vector<T> u2(v2);
    if (u1.size() != u2.size()) return false;
    std::sort(u1.begin(), u1.end());
    std::sort(u2.begin(), u2.end());
    return u1 == u2;
  }

  /**
   * @return 'true' if the two std::vectors contains the same elements, <b>in the same frequency</b>, whatever their order in the container.
   *
   * @warning The two input std::vectors will be sorted.
   *
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   */
  template<class T>
  static bool haveSameElements(std::vector<T>& v1, std::vector<T>& v2)
  {
    if (v1.size() != v2.size()) return false;
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    return v1 == v2;
  }

  /**
   * @return 'true' if a the input std::vector contains the given element.
   * @param vec The std::vector to check.
   * @param el The element to look for.
   */
  template<class T>
  static bool contains(const std::vector<T>& vec, T el)
  {
    for (size_t i = 0; i < vec.size(); i++)
    {
      if (vec[i] == el) return true;
    }
    return false;
  }

  /**
   * @return 'true' if a the first std::vector contains all elements of the second std::vector.
   *
   * @warning The two input std::vectors will be sorted.
   *
   * @param v1 The first std::vector to check.
   * @param v2 The second std::vector to check.
   */
  template<class T>
  static bool containsAll(std::vector<T>& v1, std::vector<T>& v2)
  {
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    size_t j = 0;
    for (size_t i = 0; i < v2.size(); i++)
    {
      if (i > 0 && v2[i] == v2[i - 1]) continue;
      while (j < v1.size() - 1 && v1[j] < v2[i]) j++;
      if (v1[j] != v2[i]) return false;
    }
    return true;
  }

  /**
   * @return A std::vector which is the union of two std::vectors passed as input.
   * Duplicate element will be removed.
   * @param vec1 Vector 1.
   * @param vec2 Vector 2.
   */
  template<class T>
  static std::vector<T> vectorUnion(const std::vector<T>& vec1, const std::vector<T>& vec2)
  {
    std::vector<T> unionEl = vec1;
    for (size_t j = 0; j < vec2.size(); j++)
    {
      if (!contains(unionEl, vec2[j]))
        unionEl.push_back(vec2[j]);
    }
    return unionEl;
  }

  /**
   * @return A std::vector which is the union of all std::vectors passed as input.
   * Duplicate element will be removed.
   * @param vecElementL A std::vector of std::vectors.
   */
  template<class T>
  static std::vector<T> vectorUnion(const std::vector< std::vector<T> >& vecElementL)
  {
    std::vector<T> unionEl;
    for (size_t i = 0; i < vecElementL.size(); i++)
    {
      for (size_t j = 0; j < vecElementL[i].size(); j++)
      {
        if (!contains(unionEl, vecElementL[i][j]))
          unionEl.push_back(vecElementL[i][j]);
      }
    }
    return unionEl;
  }

  /**
   * @return A std::vector which is the intersection of two std::vectors passed as input.
   * @param vec1 Vector 1.
   * @param vec2 Vector 2.
   */
  template<class T>
  static std::vector<T> vectorIntersection(const std::vector<T>& vec1, const std::vector<T>& vec2)
  {
    std::vector<T> interEl;
    for (size_t i = 0; i < vec1.size(); i++)
    {
      if (contains(vec2, vec1[i])) interEl.push_back(vec1[i]);
    }
    return interEl;
  }

  /**
   * @return A std::vector which is the intersection of all std::vectors passed as input.
   * @param vecElementL A std::vector of std::vectors.
   */
  template<class T>
  static std::vector<T> vectorIntersection(const std::vector< std::vector<T> >& vecElementL)
  {
    if (vecElementL.size() == 1) return vecElementL[0];
    std::vector<T> interEl;
    if (vecElementL.size() == 0) return interEl;
    for (size_t i = 0; i < vecElementL[0].size(); i++)
    {
      bool test = true;
      for (size_t j = 1; test && j < vecElementL.size(); j++)
      {
        if (!contains(vecElementL[j], vecElementL[0][i])) test = false;
      }
      if (test) interEl.push_back(vecElementL[0][i]);
    }
    return interEl;
  }

  /**
   * @brief Append the content of a std::vector to another one.
   * @param vec1 Vector 1.
   * @param vec2 Vector 2.
   */
  template<class T>
  static void append(std::vector<T>& vec1, const std::vector<T>& vec2)
  {
    vec1.insert(vec1.end(), vec2.begin(), vec2.end());
    // for(size_t i = 0; i < vec2.size(); i++)
    // {
    //  vec1.push_back(vec2[i]);
    // }
  }

  /**
   * @brief Prepend the content of a std::vector to another one.
   * @param vec1 Vector 1.
   * @param vec2 Vector 2.
   */
  template<class T>
  static void prepend(std::vector<T>& vec1, const std::vector<T>& vec2)
  {
    vec1.insert(vec1.begin(), vec2.begin(), vec2.end());
  }


  /**
   * @return A single std::vector made of the concatenation of the std::vectors passed as input.
   * @param vecElementL A std::vector of std::vectors.
   */
  template<class T>
  static std::vector<T> append(const std::vector< std::vector<T> >& vecElementL)
  {
    if (vecElementL.size() == 1) return vecElementL[0];
    std::vector<T> v;
    if (vecElementL.size() == 0) return v;
    for (size_t i = 0; i < vecElementL[0].size(); i++)
    {
      v.push_back(vecElementL[0][i]);
    }
    return v;
  }

  /**
   * @brief Extend the content of a std::vector with another one. Only the elements not present in the first vector will be added.
   * @param vec1 Vector 1.
   * @param vec2 Vector 2.
   */
  template<class T>
  static void extend(std::vector<T>& vec1, const std::vector<T>& vec2)
  {
    for (size_t i = 0; i < vec2.size(); i++)
    {
      if (!contains(vec1, vec2[i]))
        vec1.push_back(vec2[i]);
    }
  }

  /**
   * @return A single std::vector made of the repetion of the std::vectors passed as input.
   * @param vec A std::vector.
   * @param n the number of repetitions
   */
  template<class T>
  static std::vector<T> rep(const std::vector<T>& vec, size_t n)
  {
    if (n == 1) return vec;
    std::vector<T> v;
    if (n == 0) return v;
    v.resize(vec.size() * n);
    for (size_t i = 0; i < v.size(); i++)
    {
      v[i] = vec[i % vec.size()];
    }
    return v;
  }

  /**
   * @brief This function returns the difference of two std::vectors.
   *
   * @warning The two input std::vectors will be sorted. As a consequence, the output std::vector will be also sorted.
   *
   * @param v1 First std::vector.
   * @param v2 Second std::vector.
   * @param v3 A std::vector to be populated with all elements in v1 that are not found in v2.
   */
  template<class T>
  static void diff(std::vector<T>& v1, std::vector<T>& v2, std::vector<T>& v3)
  {
    if (v2.size() == 0) append(v3, v1);
    std::sort(v1.begin(), v1.end());
    std::sort(v2.begin(), v2.end());
    size_t j = 0;
    for (size_t i = 0; i < v1.size(); i++)
    {
      if (i > 0 && v1[i] == v1[i - 1]) continue;
      while (j < v2.size() - 1 && v2[j] < v1[i]) j++;
      if (v2[j] != v1[i]) v3.push_back(v1[i]);
    }
  }

  /**
   * @brief Test function.
   * @return true if all tests are passed.
   */
  static bool test();
};
} // end of namespace bpp.

#endif  // _VECTORTOOLS_H_

