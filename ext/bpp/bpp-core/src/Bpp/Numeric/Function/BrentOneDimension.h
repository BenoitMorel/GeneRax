//
// File: BrentOneDimension.h
// Created by: Julien Dutheil
// Created on: Mon Nov 17 11:45:58 2003
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

#ifndef _BRENTONEDIMENSION_H_
#define _BRENTONEDIMENSION_H_

#include "AbstractOptimizer.h"

namespace bpp
{

/**
 * @brief Brent's optimization for one parameter.
 *
 * A description of the algorithm can be found in:
 * <pre>
 * NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
 * (ISBN 0-521-43108-5)
 * </pre>
 * or there: <a href="http://en.wikipedia.org/wiki/Brent's_method">http://en.wikipedia.org/wiki/Brent's_method</a>.
 */
class BrentOneDimension:
  public AbstractOptimizer
{
	public:
		class BODStopCondition:
      public AbstractOptimizationStopCondition
		{
      public:
				BODStopCondition(BrentOneDimension* bod):
          AbstractOptimizationStopCondition(bod) {
            tolerance_ = bod->tol2;
            burnin_ = 3;
          }
				virtual ~BODStopCondition() {}

        BODStopCondition* clone() const { return new BODStopCondition(*this); } 
			
			public:
				bool isToleranceReached() const;
				double getCurrentTolerance() const;
		};
	
	friend class BODStopCondition;
	
	protected:
		double a, b, d, e, etemp, fu, fv, fw, fx, p, q, r, tol1, tol2, u, v, w, x, xm;
		double _xinf, _xsup;
    bool isInitialIntervalSet_;

	public:
		BrentOneDimension(Function* function = 0);
		virtual ~BrentOneDimension() {}

    BrentOneDimension* clone() const { return new BrentOneDimension(*this); }
	
	public:		
		
		/**
		 * @name The Optimizer interface.
		 *
		 * @{
		 */
		
		/**
		 * @brief Initialize optimizer.
		 *
		 * Brent's algorithm needs 2 initial guesses, so you must call the
		 * setInitialInterval() method first. This function actually performs:
     * <ul>
		 * <li>Parameter list actualisation;</li>
		 * <li>Initial bracketting;</li>
		 * <li>Function evaluation count reseting.</li>
     * </ul>
		 */
    double optimize() throw (Exception); //redefinition
		/** @} */
		
    void doInit(const ParameterList& params) throw (Exception);
		
    double doStep() throw (Exception);
	
	public:

		/**
		 * @name Specific method
		 *
		 * @{
		 */
		
		/**
		 * @brief Set intial search interval.
		 *
		 * @param inf Minimum value.
		 * @param sup Maximum value.
		 */
		void setInitialInterval(double inf, double sup);
		/** @} */

    /**
     * @return 'true' if the initial interval has been correctly set.
     */
    bool isInitialIntervalSet() const { return isInitialIntervalSet_; }
	
	public:
		
	static double ZEPS;

};

} //end of namespace bpp.

#endif	//_BRENTONEDIMENSION_H_

