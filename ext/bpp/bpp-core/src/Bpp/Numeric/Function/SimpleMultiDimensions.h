//
// File: SimpleMultiDimensions.h
// Created by: Julien Dutheil
// Created on: Tue Nov 16 17:51 2004
//

/*
Copyright or © or Copr. Bio++ Development Team, (November 19, 2004)

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

#ifndef _SIMPLEMULTIDIMENSIONS_H_
#define _SIMPLEMULTIDIMENSIONS_H_

#include "AbstractOptimizer.h"
#include "BrentOneDimension.h"

namespace bpp
{

/**
 * @brief This Optimizer is a very simple multi-dimensions optimizer, calling
 * a one dimensional optimizer on each parameter.
 *
 * The one-dimensional optimizer used is BrentOneDimension.
 * Consider using PowellMultiDimensions optimizer for a more efficient modified version of the algorithm.
 */
class SimpleMultiDimensions:
  public AbstractOptimizer
{
	private:

		size_t nbParams_;

		BrentOneDimension optimizer_; // One dimensional optimizer.

	public:

		SimpleMultiDimensions(Function* function);

		virtual ~SimpleMultiDimensions() {}

    SimpleMultiDimensions* clone() const { return new SimpleMultiDimensions(*this); }

	public:
		/**
		 * @name The Optimizer interface.
		 *
		 * @{
		 */
		void setFunction(Function * function);
		/** @} */
		
		void doInit(const ParameterList& params) throw (Exception);

    double doStep() throw (Exception);

    /**
     * @return The optimizer used to optimize each parameter.
     */
    Optimizer& getOneDimensionOptimizer() { return optimizer_; } 
    /**
     * @return The optimizer used to optimize each parameter.
     */
    const Optimizer& getOneDimensionOptimizer() const { return optimizer_; } 
};

} //end of namespace bpp.

#endif //_SIMPLEMULTIDIMENSIONS_H_

