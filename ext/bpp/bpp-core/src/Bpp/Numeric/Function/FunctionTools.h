//
// File: FunctionTools.h
// Created by: Julien Dutheil
// Created on: Mon Apr 13 10:00 2009
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

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

#ifndef _FUNCTIONTOOLS_H_
#define _FUNCTIONTOOLS_H_

#include "Functions.h"
#include "../VectorTools.h"

namespace bpp
{

/**
 * @brief This class is a data structure to specify a set of parameter values (most likely for evaluation by a Function)
 *
 * @see FunctionTools
 */
class ParameterGrid
{
  private:
    std::vector<std::string> names_;
    VVdouble grid_;

  public:
    ParameterGrid(): names_(), grid_() {}
    virtual ~ParameterGrid() {}

  public:
    /**
     * @brief Add a new dimension (parameter name + corresponding values).
     *
     * @param name The name of the dimension (parameter name).
     * @param values The values the parameter will take.
     * @throw Exception in case the dimension is note valid (duplicated parameter name for instance).
     */
    void addDimension(const std::string& name, const Vdouble& values) throw (Exception);

    const std::vector<std::string>& getDimensionNames() const { return names_; }
    
    const std::string& getDimensionName(unsigned int i) const throw (IndexOutOfBoundsException)
    {
      if (i >= names_.size()) throw IndexOutOfBoundsException("ParameterGrid::getDimensionName().", i, 0, names_.size()-1);
      return names_[i];
    }

    size_t getNumberOfDimensions() const { return names_.size(); }
    
    /**
     * @return The total number of points in the grid, that is the product of all dimension sizes.
     */
    size_t getTotalNumberOfPoints() const;

    const VVdouble& getPoints() const { return grid_; }
    const Vdouble& getPointsForDimension(unsigned int i) const throw (IndexOutOfBoundsException);
    const Vdouble& getPointsForDimension(const std::string& name) const throw (Exception);
};

/**
 * @brief This class contains static methods to deal with Function objects.
 */
class FunctionTools
{
  public:
    /**
     * @brief Evaluates a function on all points in a given grid.
     *
     * @param function The function to use for the evaluation.
     * @param grid     The grid defining the set of points to evaluate.
     * @return A pointer toward a dynamically created vector of vector
     * of doubles. Each row correpsonds to a combination of parameters
     * and the corresponding function value. There is hence one column
     * per parameter, and one additional column containing the
     * corresponding function evaluations. When DataTable supports
     * different column type, we will probably return a DataTable instead.
     * @throw Exception If the parameter names in the grid do not match
     * the ones in the function, or a constraint is matched, etc.
     */
    static VVdouble* computeGrid(
        Function& function,
        const ParameterGrid& grid) throw (Exception);
};

} //end of namespace bpp

#endif //_FUNCTIONTOOLS_H_

