//
// File: MetaOptimizer.h
// Created by: Julien Dutheil
// Created on: Fri Oct 12 16:05 2007
// From file: NewtonBrentMetaOptimizer.h
// Created on: Tue Nov 17 17:22 2004
// 
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

#ifndef METAOPTIMIZER_H__
#define METAOPTIMIZER_H__

#include "AbstractOptimizer.h"

// From the STL:
#include <vector>

namespace bpp
{

  /**
   * @brief Provide a list of optimizer and corresponding options to be used with the MetaOptimizer class.
   */
  class MetaOptimizerInfos:
    public virtual Clonable
  {
  public:
    static std::string IT_TYPE_STEP;
    static std::string IT_TYPE_FULL;

  private:
    std::vector<std::string> names_;
    std::vector<Optimizer*> optimizers_;
    std::vector< std::vector<std::string> > parameterNames_;
    std::vector<unsigned short> derivatives_;
    std::vector<std::string> itTypes_;

  public:
    MetaOptimizerInfos() : names_(), optimizers_(), parameterNames_(), derivatives_(), itTypes_() {}
    MetaOptimizerInfos(const MetaOptimizerInfos& infos) :
      names_(infos.names_),
      optimizers_(infos.optimizers_),
      parameterNames_(infos.parameterNames_),
      derivatives_(infos.derivatives_),
      itTypes_(infos.itTypes_)
    {
      for (unsigned int i = 0; i < optimizers_.size(); i++)
        optimizers_[i] = dynamic_cast<Optimizer*>(infos.optimizers_[i]->clone());
    }

    MetaOptimizerInfos& operator=(const MetaOptimizerInfos& infos)
    {
      names_          = infos.names_;
      optimizers_     = infos.optimizers_;
      parameterNames_ = infos.parameterNames_;
      derivatives_    = infos.derivatives_;
      itTypes_        = infos.itTypes_;
      for (unsigned int i = 0; i < optimizers_.size(); i++)
        optimizers_[i] = dynamic_cast<Optimizer *>(infos.optimizers_[i]->clone());
      return *this;
    }

    virtual ~MetaOptimizerInfos()
    {
      for(unsigned int i = 0; i < optimizers_.size(); i++)
        delete optimizers_[i];
    }

  public:
    MetaOptimizerInfos* clone() const { return new MetaOptimizerInfos(*this); }

  public:
    /**
     * @brief Add a new optimizer to the set.
     *
     * @param name the name of the optimizer. It is used for display only.
     * @param optimizer A pointer toward the optimizer to add. The set will own the underlying object, which will be destroyed together with the set.
     * @param params A list of parameter names to optimize with this optimizer.
     * @param derivatives 0, 1 or 2: does this parameter use no, first order or second order derivatives?
     * @param type For each optimization step, shall we perform a full optimization with this optimizer or only one step?
     */
    virtual void addOptimizer(const std::string & name, Optimizer * optimizer, const std::vector<std::string> & params, unsigned short derivatives = 0, const std::string & type = IT_TYPE_STEP)
    {
      names_.push_back(name);
      optimizers_.push_back(optimizer);
      parameterNames_.push_back(params);
      derivatives_.push_back(derivatives);
      itTypes_.push_back(type);
    }

    /**
     * @return The display name of the ith optimizer in the set.
     */
    virtual const std::string& getName(size_t i) const { return names_[i]; }
    
    /**
     * @return The ith optimizer in the set.
     */
    virtual Optimizer* getOptimizer(size_t i) { return optimizers_[i]; }
    /**
     * @return The ith optimizer in the set.
     */
    virtual const Optimizer* getOptimizer(size_t i) const { return optimizers_[i]; }

    /**
     * @return The parameter names associated to the ith optimizer in the set.
     */
    virtual std::vector<std::string>& getParameterNames(size_t i) { return parameterNames_[i]; }
    /**
     * @return The parameter names associated to the ith optimizer in the set.
     */
    virtual const std::vector<std::string>& getParameterNames(size_t i) const { return parameterNames_[i]; }

    /**
     * @return The type of iteration to perform for the ith optimizer in the set.
     */
    virtual std::string& getIterationType(size_t i) { return itTypes_[i]; }
    /**
     * @return The type of iteration to perform for the ith optimizer in the set.
     */
    virtual const std::string& getIterationType(size_t i) const { return itTypes_[i]; }

    /**
     * @return True if the ith optimizer in the set requires first order derivatives.
     */
    virtual bool requiresFirstOrderDerivatives(size_t i) const { return derivatives_[i] > 0; }  
    /**
     * @return True if the ith optimizer in the set requires second order derivatives.
     */
    virtual bool requiresSecondOrderDerivatives(size_t i) const { return derivatives_[i] > 1; }  

    /**
     * @return The number of optimizers in the set.
     */
    virtual size_t getNumberOfOptimizers() const { return optimizers_.size(); }
  };

  /**
   * @brief Meta-optimizer.
   *
   * This optimizer uses a set of optimizers to applyied sequentially on distinct parameters.
   * The set of optimizers is fully specified by a MetaOptimizerInfos object.
   * 
   * To decrease the optimization time, the precision of the optimizers can be increased progressively:
   * if @f$\varepsilon@f$ is the final precision required, one may consider using a precision increment of @f$\sigma=\log_10(\varepsilon/n)@f$, where @f$n@f$ is the number of progressive steps.
   * During the first step optimization step, the precisions of type 1 and 2 optimizers are set to @f$10^{\sigma}@f$, @f$10^{2\sigma}@f$ for step 2, ... until precision @f$10^{n\sigma}=\varepsilon@f$ at step @f$n@f$ and later.
   * This saves some time spending in the first steps of the estimation.
   * The number of steps @f$n@f$ is set in the constructor of the optimizer.
   *
   * This optimizer can be used with numerical derivatives.
   * 
   * @see MetaOptimizerInfos.
   */
  class MetaOptimizer:
    public AbstractOptimizer
  {
  private:
    MetaOptimizerInfos* optDesc_;
    std::vector<ParameterList> optParameters_;
    std::vector<size_t> nbParameters_;
    unsigned int n_;
    double precisionStep_;
    unsigned int stepCount_;
    double initialValue_;
		
  public:
    /**
     * @brief Build a new MetaOptimizer object.
     *
     * @param function The function to be optimized.
     * @param desc     A MetaOptimizerInfos object describing the optimizers to use.
     *                 The optimizer will own the instance of the MetaOptimizerInfos object.
     * @param n        The number of progressive steps to use in optimization).
     */
    MetaOptimizer(Function* function, MetaOptimizerInfos* desc, unsigned int n = 1);

    virtual ~MetaOptimizer();

    MetaOptimizer(const MetaOptimizer& opt);

    MetaOptimizer& operator=(const MetaOptimizer& opt);

    MetaOptimizer* clone() const { return new MetaOptimizer(*this); }

  public:
		
    void setFunction(Function* function)
    { 
      AbstractOptimizer::setFunction(function);
      for(unsigned int i = 0; i < optDesc_->getNumberOfOptimizers(); i++)
        optDesc_->getOptimizer(i)->setFunction(function);
    }

    void doInit(const ParameterList& parameters) throw (Exception);
    
    double doStep() throw (Exception);

    /**
     * @return The MetaOptimizerInfos object associated to this optimizer.
     */
    MetaOptimizerInfos* getOptimizers() { return optDesc_; }
    
    /**
     * @return The MetaOptimizerInfos object associated to this optimizer.
     */
    const MetaOptimizerInfos* getOptimizers() const { return optDesc_; }

  };

} //end of namespace bpp.

#endif //METAOPTIMIZER_H__

