//
// File: ConjugateGradientMultiDimensions.h
// Created by: Julien Dutheil
// Created on: Wed Apr 11 16:51 2007
//

/*
  Copyright or © or Copr. CNRS, (November 19, 2004)

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

#ifndef _CONJUGATEGRADIENTMULTIDIMENSIONS_H_
#define _CONJUGATEGRADIENTMULTIDIMENSIONS_H_

#include "AbstractOptimizer.h"
#include "BrentOneDimension.h"
#include "DirectionFunction.h"

namespace bpp
{

  /**
   * @brief Conjugate gradient optimization method.
   *
   * A description of the algorithm can be found in:
   * <pre>
   * NUMERICAL RECIPES IN C: THE ART OF SCIENTIFIC COMPUTING
   * (ISBN 0-521-43108-5)
   * </pre>
   * or there:
   * <a href="http://en.wikipedia.org/wiki/Conjugate_gradient">http://en.wikipedia.org/wiki/Conjugate_gradient</a>.
   */
  class ConjugateGradientMultiDimensions:
    public AbstractOptimizer
  {
  protected:
    BrentOneDimension optimizer_; //One dimensional optimizer.
    std::vector<double> xi_, h_, g_;
    DirectionFunction f1dim_;

  public:

    ConjugateGradientMultiDimensions(DerivableFirstOrder* function);

    virtual ~ConjugateGradientMultiDimensions() {}

    ConjugateGradientMultiDimensions* clone() const { return new ConjugateGradientMultiDimensions(*this); }

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
    
  };

} //end of namespace bpp.

#endif //_CONJUGATEGRADIENTMULTIDIMENSIONS_H_

