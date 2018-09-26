//
// File: OptimizationStopCondition.h
// Created by: Julien Dutheil
// Created on: Tue Dec 23 11:51:31 2003
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

#ifndef _OPTIMIZATIONSTOPCONDITION_H_
#define _OPTIMIZATIONSTOPCONDITION_H_

#include "../ParameterList.h"

using namespace std;

namespace bpp
{

  class Optimizer;
	
  /******************************************************************************/

  /**
   * @brief Interface for otimization stop condition objet.
   *
   * Classes implementing the OptimizationStopCondition interface
   * provides the isToleranceReached() function that tells when
   * the optimization process reached a given tolerance parameter.
   * This parameter may be set or retrieve using the setTolerance()
   * and getTolerance() functions.
   *
   * OptimizationStopCondition implementations may be general and work on
   * parameter (@see ParametersStopCondition) or function (@see FunctionStopCondition) values,
   * or be specific to a given optimization method.
   */
  class OptimizationStopCondition:
    public virtual Clonable
  {
  public:
    OptimizationStopCondition() {}
    virtual ~OptimizationStopCondition() {}

    OptimizationStopCondition * clone() const = 0;
	
  public:

    /**
     * @return The optimizer to which this instance belongs to. 
     */
    virtual const Optimizer* getOptimizer() const = 0;
    /**
     * @brief Set the optimizer attached to this instance.
     *
     * @param optimizer The optimizer to which this instance belongs to. 
     */
    virtual void setOptimizer(const Optimizer* optimizer) = 0;

    /**
     * @brief Initialize the condition.
     */
    virtual void init() = 0;

    /**
     * @brief Tell if the we reached the desired tolerance with a given 
     * new set of estimates.
     *
     * The new parameter list is compared to the last estimates,
     * and the lastParameterEstimates list is actulaized with the newParameters list.
     *
     * @return True if the tolerance level is reached.
     */
    virtual bool isToleranceReached() const = 0;
		
    /**
     * @brief Set the tolerance parameter.
     *
     * @param tolerance The tolerance parameter.
     */	
    virtual void setTolerance(double tolerance) = 0;

    /**
     * @brief Get the tolerance parameter.
     *
     * @return The tolerance parameter.
     */	
    virtual double getTolerance() const = 0;
    
    /**
     * @brief Get the current tolerance.
     *
     * This is computed from the last check performed.
     * Initially, it is equal to the tolerance parameter.
     *
     * @return The current tolerance achieved.
     */	
    virtual double getCurrentTolerance() const = 0;
  };

  /******************************************************************************/

  /**
   * @brief Partial implementation of the OptimizationStopCondition interface.
   *
   * This class provides:
   * - A pointer toward the Optimizer this objet deals with,
   * - A tolerance value,
   * - A counter of the number of calls toward the isToleranceReached() function,
   * - A burnin function, that prohibe the optimization to stop prematurely.
   */   
  class AbstractOptimizationStopCondition:
    public virtual OptimizationStopCondition
  {
  protected:
    const Optimizer* optimizer_;
    double tolerance_;

    /**
     * @brief Count the number of times the isToleranceReached() function
     * has been called.
     */
    mutable double callCount_;
	
    int burnin_;
	
  public:
    AbstractOptimizationStopCondition(const Optimizer* optimizer):
        optimizer_(optimizer),
        tolerance_(0.000001),
        callCount_(0),
        burnin_(0) {}
    
    AbstractOptimizationStopCondition(const Optimizer* optimizer, double tolerance):
        optimizer_(optimizer),
        tolerance_(tolerance),
        callCount_(0),
        burnin_(0) {}

    AbstractOptimizationStopCondition(const Optimizer* optimizer, int burnin):
        optimizer_(optimizer),
        tolerance_(0.000001),
        callCount_(0),
        burnin_(burnin) {}

    AbstractOptimizationStopCondition(const Optimizer* optimizer, double tolerance, int burnin):
        optimizer_(optimizer),
        tolerance_(tolerance),
        callCount_(0),
        burnin_(burnin) {}

    AbstractOptimizationStopCondition(const AbstractOptimizationStopCondition& aosc):
        optimizer_(aosc.optimizer_),
        tolerance_(aosc.tolerance_),
        callCount_(aosc.callCount_),
        burnin_(aosc.burnin_) {}
	
    AbstractOptimizationStopCondition& operator=(const AbstractOptimizationStopCondition& aosc)
    {
      optimizer_ = aosc.optimizer_;
      tolerance_ = aosc.tolerance_;
      callCount_ = aosc.callCount_;
      burnin_    = aosc.burnin_;
      return *this;
    }

    virtual ~AbstractOptimizationStopCondition() {}

  public:
    const Optimizer* getOptimizer() const { return optimizer_; }
    void setOptimizer(const Optimizer* optimizer) { optimizer_ = optimizer; }
    void setTolerance(double tolerance) { tolerance_ = tolerance; }
    double getTolerance() const { return tolerance_; }
    void init() { resetCounter(); }
    virtual void resetCounter() { callCount_ = 0; }
    virtual void setBurnin(int burnin) { burnin_ = burnin; }
    virtual int getBurnin() const { return burnin_; }

  };
	
  /******************************************************************************/

  /**
   * @brief Stop condition on parameters.
   *
   * This stops the optimization when \f$\forall i\; |\lambda_{i,t}-\lambda_{i,t-1}| \leq \mbox{tolerance}\f$,
   * where \f$\lambda_{i, t}\f$ is the value of the ith parameter at iteration \f$t\f$,
   * and \f$\lambda_{i, t-1}\f$ is the value of the ith parameter at iteration \f$t-1\f$.
   */
  class ParametersStopCondition:
    public AbstractOptimizationStopCondition
  {
  private:

    /**
     * @brief The last estimates of the parameters.
     *
     * This is used by the isToleranceReached() method.
     */
    mutable ParameterList lastParametersEstimates_;
		
    /**
     * @brief The new estimates of the parameters.
     *
     * This is used by the isToleranceReached() method.
     */
    mutable ParameterList newParametersEstimates_;
	
  public:
    ParametersStopCondition(const Optimizer* optimizer);
    ParametersStopCondition(const Optimizer* optimizer, double tolerance);
    ParametersStopCondition(const Optimizer* optimizer, int burnin);
    ParametersStopCondition(const Optimizer* optimizer, double tolerance, int burnin);
		
    virtual ~ParametersStopCondition() {}

    ParametersStopCondition* clone() const { return new ParametersStopCondition(*this); }
	
  public:
    void init();

    bool isToleranceReached() const;
    
    double getCurrentTolerance() const;
  };

  /******************************************************************************/

  /**
   * @brief Stop condition on function value.
   *
   * This stops the optimization when \f$|f\left(\left\{\lambda_{i,t}\right\}\right)-f\left(\left\{\lambda_{i,t-1}\right\}\right)| \leq \mbox{tolerance}\f$,
   * where \f$f\left(\left\{\lambda_{i, t}\right\}\right)\f$ is the value of the function given the parameter values at iteration \f$t\f$,
   * and \f$f\left(\left\{\lambda_{i, t-1}\right\}\right)\f$ is the value of the function given the parameter velues at iteration \f$t-1\f$.
   */
  class FunctionStopCondition:
    public AbstractOptimizationStopCondition
  {
  private:
    /**
     * @brief The last value of the function.
     *
     * This is used by the isToleranceReached() method.
     */
    mutable double lastFunctionValue_;
		
    /**
     * @brief The new value of the function.
     *
     * This is used by the isToleranceReached() method.
     */
    mutable double newFunctionValue_;
	
  public:
    FunctionStopCondition(const Optimizer* optimizer);
    FunctionStopCondition(const Optimizer* optimizer, double tolerance);
    FunctionStopCondition(const Optimizer* optimizer, int burnin);
    FunctionStopCondition(const Optimizer* optimizer, double tolerance, int burnin);
		
    virtual ~FunctionStopCondition();

    FunctionStopCondition* clone() const { return new FunctionStopCondition(*this); }
	
  public:
    void init();
    bool isToleranceReached() const;
    double getCurrentTolerance() const;

  };

} //end of namespace bpp.

#endif	//_OPTIMIZATIONSTOPCONDITION_H_

