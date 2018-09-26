//
// File: FunctionTools.cpp
// Created by: Julien Dutheil
// Created on: Mon Apr 13 10:47 2009
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

#include "FunctionTools.h"
#include "../../App/ApplicationTools.h"

using namespace bpp;

//From the STL;
#include <algorithm>
using namespace std;

void ParameterGrid::addDimension(const std::string& name, const Vdouble& values) throw (Exception)
{
  if (find(names_.begin(), names_.end(), name) != names_.end()) throw Exception("ParameterGrid::addDimension(). A dimension with name '" + name + "' already exists in the grid.");
  if (values.size() == 0) throw Exception("ParameterGrid::addDimension(). Empty vector given! The dimension should at least contain one point.");
  names_.push_back(name);
  grid_.push_back(values);
}
 
const Vdouble& ParameterGrid::getPointsForDimension(const std::string& name) const throw (Exception)
{
  for(unsigned int i = 0; i < names_.size(); i++)
    if (names_[i] == name)
      return grid_[i];
  throw Exception("ParameterGrid::getPointsForDimension(). No dimension with name '" + name + "' was found in the grid.");
}

const Vdouble& ParameterGrid::getPointsForDimension(unsigned int i) const throw (IndexOutOfBoundsException)
{
  if (i >= names_.size()) throw IndexOutOfBoundsException("ParameterGrid::getPointsForDimension().", i, 0, names_.size() - 1);
  return grid_[i];
}

size_t ParameterGrid::getTotalNumberOfPoints() const
{
  if (grid_.size() == 0) return 0;
  size_t n = 1;
  for (size_t i = 0; i < grid_.size(); i++)
    n *= grid_[i].size();
  return n;
}

VVdouble* FunctionTools::computeGrid(
    Function& function,
    const ParameterGrid& grid) throw (Exception)
{
  //Init stuff...
  size_t n = grid.getNumberOfDimensions();
  VVdouble* data = new VVdouble();
  if(n == 0) return data; //Empty data table returned.

  VVdouble points = grid.getPoints();

  //Get the parameter list. this may throw an exception if the grid does not
  //match the function parameters...
  ParameterList pl = function.getParameters().subList(grid.getDimensionNames());
  for(unsigned int i = 0; i < n; i++)
    pl.setParameterValue(grid.getDimensionName(i), grid.getPointsForDimension(i)[0]);

  //Iterate over all dimensions:
  unsigned int currentDimension = 0;
  vector<unsigned int> currentPointInDimension(n);
  vector<double> row(n + 1);
  size_t nbPoints = grid.getTotalNumberOfPoints();
  ApplicationTools::displayMessage("Computing likelihood profile...");
  for (unsigned int i = 0; true ; i++)
  {
    ApplicationTools::displayGauge(i, nbPoints - 1, '=');
    //We start by adding the current point to the table:
    for (unsigned int j = 0; j < n; j++)
      row[j] = pl[j].getValue();
    row[n] = function.f(pl);
    data->push_back(row);

    //Now increment iterator:
    bool dimensionChanged = false;
    while (currentDimension < n && currentPointInDimension[currentDimension] == points[currentDimension].size() - 1)
    {
      currentDimension++;
      dimensionChanged = true;
    }
    //Stopping condition:
    if (currentDimension == n) break;

    currentPointInDimension[currentDimension]++;
    if (dimensionChanged)
    {
      for (unsigned int j = 0; j < currentDimension; j++)
        currentPointInDimension[j] = 0;
      currentDimension = 0;
    }
   
    //Set the new parameter value:
    for (unsigned int j = 0; j < points.size(); j++)
    {
      pl.setParameterValue(grid.getDimensionName(j), points[j][currentPointInDimension[j]]);
    }
  }
  ApplicationTools::displayMessage("\n");
  //and we are done:
  return data;
}

