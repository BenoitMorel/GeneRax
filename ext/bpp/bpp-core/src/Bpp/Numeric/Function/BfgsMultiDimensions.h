//
// File: BfgsMultiDimensions.h
// Created by: Laurent Guéguen
// Created on: Dec 16 13:49 2010
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

#ifndef _BFGSMULTIDIMENSIONS_H_
#define _BFGSMULTIDIMENSIONS_H_

#include "AbstractOptimizer.h"
#include "DirectionFunction.h"
#include "../VectorTools.h"

namespace bpp
{

  /**
   * @brief Broyden–Fletcher–Goldfarb–Shanno (BFGS) optimization method.
   *
   * with a modification on the bounds taken from:
   *  An active set limited memory BFGS algorithm for large-scale bound
   *    constrained optimization, Yunhai Xiao & Dong-Hui Li. Math Meth
   *    Oper Res (2008) 67:443­454
   */
  
  class BfgsMultiDimensions:
    public AbstractOptimizer
  {
  protected:
    //double gtol_;
    double slope_;

    // vectors of the Lower & Upper bounds of the parameters
    Vdouble Up_, Lo_;

    mutable Vdouble p_, gradient_, xi_, dg_, hdg_;
    mutable VVdouble hessian_;

    mutable DirectionFunction f1dim_;

    
  public:

    BfgsMultiDimensions(DerivableFirstOrder* function);

    virtual ~BfgsMultiDimensions() {}

    BfgsMultiDimensions* clone() const { return new BfgsMultiDimensions(*this); }

  public:

    /**
     * @name From AbstractOptimizer.
     *
     * @{
     */
    const DerivableFirstOrder* getFunction() const
    {
      return dynamic_cast<const DerivableFirstOrder*>(AbstractOptimizer::getFunction());
    }
    DerivableFirstOrder* getFunction()
    {
      return dynamic_cast<DerivableFirstOrder*>(AbstractOptimizer::getFunction());
    }
    void doInit(const ParameterList& params) throw (Exception);

    double doStep() throw (Exception);
    /** @} */

    void getGradient(std::vector<double>& gradient) const;

  protected:
    DerivableFirstOrder* getFunction_()
    {
      return dynamic_cast<DerivableFirstOrder*>(AbstractOptimizer::getFunction_());
    }

  private:
    /**
     * Sets the direction for linesearch in case of bounds
     * To be used after gradient_ & pi_ are computed
     */ 
    void setDirection();
    
  };

} //end of namespace bpp.

#endif //_BFGSMULTIDIMENSIONS_H_

