//
// File: NewtonBacktrackOneDimension.h
// Created by: Laurent Guéguen
// Created on: jeudi 16 décembre 2010, à 15h 43
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

#ifndef _NEWTONBACKTRACKONEDIMENSION_H_
#define _NEWTONBACKTRACKONEDIMENSION_H_

#include "AbstractOptimizer.h"

namespace bpp
{

  /**
   * @brief Newton's backtrack nearly optimization for one parameter.
   *
   *
   * This algorithm looks for a value that 'sufficiently low' enough
   * near the minimum of the function, but does not look for the
   * minimum. It needs the first derivative of a function.
   *
   *  Search a 'sufficiently low' value for a function in a given
   *  direction.
   *
   * This function implements the algorithm described for example in page 385 of
   *
   * <pre>
   * NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
   * (ISBN 0-521-43108-5)
   * </pre>
   */

  class NewtonBacktrackOneDimension:
    public AbstractOptimizer
  {
  public:
    class NBODStopCondition:
      public AbstractOptimizationStopCondition
    {
    public:
      NBODStopCondition(NewtonBacktrackOneDimension* bod):
        AbstractOptimizationStopCondition(bod) {}
      virtual ~NBODStopCondition() {}
      
      NBODStopCondition* clone() const { return new NBODStopCondition(*this); } 
      
    public:
      void init() {}
      bool isToleranceReached() const { return false; }
      double getCurrentTolerance() const { return 0; }
    };
    
    friend class NBODStopCondition;
    
  public:

    /*
     *@brief Constructor
     @param function  a pointer to the function
     @param slope the slope of the backtrack
     @param test the inverse factor on the stop tolerance
     *
     */
    
    NewtonBacktrackOneDimension(Function* function, double slope, double test);
    virtual ~NewtonBacktrackOneDimension() {}

    NewtonBacktrackOneDimension* clone() const { return new NewtonBacktrackOneDimension(*this); }

  private:
    double fold_, f_, a_, alam_, alamin_, alam2_, b_, disc_, f2_, rhs1_, rhs2_, slope_, test_, tmplam_;
    
  public:
		
    void doInit(const ParameterList& params) throw (Exception);
		
    double doStep() throw (Exception);

  protected:
    DerivableFirstOrder* getFunction_()
    {
      return dynamic_cast<DerivableFirstOrder*>(AbstractOptimizer::getFunction_());
    }
    
  };

} //end of namespace bpp.

#endif	//_NEWTONBACKTRACKONEDIMENSION_H_

