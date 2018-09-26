//
// File: DownhillSimplexMethod.h
// Created by: Julien Dutheil
// Created on: Tue Nov  4 17:10:05 2003
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

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

#ifndef _DOWNHILLSIMPLEXMETHOD_H_
#define _DOWNHILLSIMPLEXMETHOD_H_

#include "AbstractOptimizer.h"
#include "../VectorTools.h"

// From the STL:
#include <cmath>

namespace bpp
{

/**
 * @brief This implements the Downhill Simplex method in multidimensions.
 *
 * A description of the algorithm can be found in:
 * <pre>
 * NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
 * (ISBN 0-521-43108-5)
 * </pre>
 * or there:
 * <a href="http://en.wikipedia.org/wiki/Nelder-Mead_method">http://en.wikipedia.org/wiki/Nelder-Mead_method</a>.
 */
class DownhillSimplexMethod:
  public AbstractOptimizer
{
  public:
    class DSMStopCondition:
      public AbstractOptimizationStopCondition
    {
      public:
        DSMStopCondition(DownhillSimplexMethod * dsm):
          AbstractOptimizationStopCondition(dsm) {}
        virtual ~DSMStopCondition() {}

        DSMStopCondition* clone() const { return new DSMStopCondition(*this); }
      
      public:
        bool isToleranceReached() const { return (getCurrentTolerance() < tolerance_); }
        double getCurrentTolerance() const;
    };
  
  friend class DSMStopCondition;
  
  private:
    class Simplex
    {
      private:
        std::vector<ParameterList> parameters_;

      public: // Class constructor and destructor:
        Simplex(): parameters_() {}
        virtual ~Simplex() {}
      
      public: // Methods:
        const ParameterList& operator[](size_t i) const { return parameters_[i]; }
        ParameterList& operator[](size_t i) { return parameters_[i]; }
        void resize(size_t size) { parameters_.resize(size); }
        size_t getDimension() const { return parameters_[0].size(); }
    };
    
  protected:
    Simplex simplex_;
    Vdouble y_;
    ParameterList pSum_;
    unsigned int iHighest_, iNextHighest_, iLowest_;
  
  public:

    /**
     * @brief Build a new Downhill Simplex optimizer.
     *
     * @param function A pointer toward an object implementing the Optimizable interface.
     */
    DownhillSimplexMethod(Function * function);
  
    virtual ~DownhillSimplexMethod() {}

    DownhillSimplexMethod* clone() const { return new DownhillSimplexMethod(*this); }
  
  public:    
    /**
     * @name The Optimizer interface.
     *
     * @{
     */
    
    /**
     * @brief Multidimensional minimization of the function function_ by the
     * downhill simplex method of Nelder and Mead.
     */
    double optimize() throw (Exception);
    /** @} */

    void doInit(const ParameterList& params) throw (Exception);
    
    double doStep() throw (Exception);
  
  protected:
    
    /**
     * @name Specific inner methods
     *
     * @{
     */
    
    /**
     * @brief Update the pSum_ variable.
     */
    ParameterList getPSum();
  
    /**
     * @brief Extrapolates by a factor fac through the face of the simplex from the high point.
     * Try the new point and replaces the high point if it is better.
     *
     * @param fac Extrapolation factor.
     * @return The value of the function for the new point.
     */
    double tryExtrapolation(double fac);

    /** @} */
};

} //end of namespace bpp.

#endif  //_DOWNHILLSIMPLEXMETHOD_H_

