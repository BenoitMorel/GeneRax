//
// File: GoldenSectionSearch.h
// Created by: Julien Dutheil
// Created on: Mon Nov 10 10:42:17 2003
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

#ifndef _GOLDENSECTIONSEARCH_H_
#define _GOLDENSECTIONSEARCH_H_

#include "AbstractOptimizer.h"

namespace bpp
{

/**
 * @brief Golden Section Search optimization algorithm for one parameter.
 *
 * A description of the algorithm can be found in:
 * <pre>
 * NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
 * (ISBN 0-521-43108-5)
 * </pre>
 * or there:
 * <a href="http://en.wikipedia.org/wiki/Golden_section_search">http://en.wikipedia.org/wiki/Golden_section_search</a>.
 */
class GoldenSectionSearch:
  public AbstractOptimizer
{	
	public:
		class GSSStopCondition:
      public AbstractOptimizationStopCondition
		{
			public:
				GSSStopCondition(GoldenSectionSearch* gss):
          AbstractOptimizationStopCondition(gss) {}
				virtual ~GSSStopCondition() {}

        GSSStopCondition* clone() const { return new GSSStopCondition(*this); }
			
			public:
				bool isToleranceReached() const;
        double getCurrentTolerance() const;
		};
	
	friend class GSSStopCondition;

	private:
		double f1, f2, x0, x1, x2, x3;
		double xinf_, xsup_;
    bool isInitialIntervalSet_;
	
	public:
		
		GoldenSectionSearch(Function* function);
		virtual ~GoldenSectionSearch() {}

    GoldenSectionSearch* clone() const { return new GoldenSectionSearch(*this); }
	
	public:
		
		/**
		 * @name The Optimizer interface.
		 *
		 * @{
		 */
		
		/**
		 * @brief Initialize optimizer.
		 *
		 * The golden section search needs 2 initial guesses, so you must call the
		 * setInitialInterval() method first. This function actually performs:
		 * <ul>
		 * <li>Parameter list actualisation;</li>
		 * <li>Initial bracketting;</li>
		 * <li>Function evaluation count reseting.</li>
		 * </ul>
		 */
		double getFunctionValue() const throw (NullPointerException);
		/** @} */
		
    void doInit(const ParameterList & params) throw (Exception);

		double doStep() throw (Exception);
	
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

  protected:
		
};

} //end of namespace bpp.

#endif	//_GOLDENSECTIONSEARCH_H_

