//
// File: NewtonOneDimension.h
// Created by: Julien Dutheil
// Created on: Thu Apr 26 14:16 2007
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

#ifndef _NEWTONONEDIMENSION_H_
#define _NEWTONONEDIMENSION_H_

#include "AbstractOptimizer.h"

namespace bpp
{

/**
 * @brief Newton's optimization for one parameter.
 */
class NewtonOneDimension:
  public AbstractOptimizer
{
	private:
    std::string _param;
    unsigned int _maxCorrection;

	public:
		NewtonOneDimension(DerivableSecondOrder* function = 0);
		virtual ~NewtonOneDimension() {}

    NewtonOneDimension* clone() const { return new NewtonOneDimension(*this); }
	
	public:		
		
    const DerivableSecondOrder* getFunction() const
    {
      return dynamic_cast<const DerivableSecondOrder*>(AbstractOptimizer::getFunction());
    }
    DerivableSecondOrder* getFunction()
    {
      return dynamic_cast<DerivableSecondOrder*>(AbstractOptimizer::getFunction());
    }

    void doInit(const ParameterList& params) throw (Exception);
		
    double doStep() throw (Exception);

    void setMaximumNumberOfCorrections(unsigned int mx) { _maxCorrection = mx; }
	
  protected:
    DerivableSecondOrder* getFunction_()
    {
      return dynamic_cast<DerivableSecondOrder*>(AbstractOptimizer::getFunction_());
    }
    
};

} //end of namespace bpp.

#endif	//_NEWTONONEDIMENSION_H_

