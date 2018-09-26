//
// File: VectorTools.cpp
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

// From Utils:
#include "../Text/TextTools.h"

#include "VectorTools.h"
using namespace bpp;

// From the STL:
#include <cmath>
#include <iostream>
using namespace std;

/******************************************************************************/

vector<double> VectorTools::breaks(const vector<double>& v, unsigned int n)
{
  vector<double> out;
  vector<double> r = VectorTools::range(v);
  double part = (r[1] - r[0]) / n;
  for (unsigned int i = 0; i < n; ++i)
  {
    out.push_back(r[0] + (part * i));
  }
  out.push_back(r[1]);
  return out;
}

/******************************************************************************/

bool VectorTools::test()
{
  vector<double> x1(5);
  vector<double> x2(5);
  x1[0] = -3.4;
  x1[1] =  1.8;
  x1[2] = -2.1;
  x1[3] = -2.5;
  x1[4] =  1.0;

  x2[0] = -5.3;
  x2[1] = -4.8;
  x2[2] =  2.7;
  x2[3] =  7.2;
  x2[4] =  0.4;

  print(x1);
  print(x2);
  double m1 = mean<double, double>(x1);
  double m2 = mean<double, double>(x2);
  double v1 = var<double, double>(x1);
  double v2 = var<double, double>(x2);
  cout << "Mean x1 = " << m1 << "\tVar x1 = " << v1 << endl;
  cout << "Mean x2 = " << m2 << "\tVar x2 = " << v2 << endl;
  cov<double, double>(x1, x2);
  cor<double, double>(x1, x2);
  cos<double, double>(x1, x2);
  shannon<double, double>(x1);
  return m1 == -0.2 && m2 == 0.04 && v1 == 6.565 && v2 == 27.603;
}

/******************************************************************************/


