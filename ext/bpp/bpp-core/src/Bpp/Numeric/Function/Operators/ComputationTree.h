//
// File: ComputationTree.h
// Created by: Laurent Guéguen
// Created on: mardi 6 décembre 2016, à 00h 07
//

/*
  Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _COMPUTATION_TREE_H_
#define _COMPUTATION_TREE_H_

#include "../../../Graph/AssociationTreeGraphImplObserver.h"

#include "Operator.h"
#include "../Functions.h"
#include <memory>

namespace bpp
{
/**
 * @brief Defines a Computation Tree based on Operators.
 *
 */  
  
  class ComputationTree:
    public AssociationTreeGlobalGraphObserver<Operator,short>
  {
  private:
    std::shared_ptr<Operator> readFormula_(const std::string& formula,const std::map<std::string, Function*>& functionNames);
    
  public:
    /*
     * @brief Tree for numerical computation given a formula (such as
     * 2*f+g), and given a map from function name (f) to real
     * function.
     * 
     */
    
    ComputationTree(const std::string& formula, const std::map<std::string, Function*>& functionNames);
    
    ComputationTree* clone() const
    {
      return new ComputationTree(*this);
    }

    double getValue() const
    {
      return getRoot()->getValue();
    }

    double getFirstOrderDerivative(const std::string& variable) const
    {
      return getRoot()->getFirstOrderDerivative(variable);
    }

    double getSecondOrderDerivative(const std::string& variable) const
    {
      return getRoot()->getSecondOrderDerivative(variable);
    }

    void readFormula(const std::string& formula,const std::map<std::string, Function*>& functionNames)
    {
      readFormula_(formula, functionNames);
    }
    
    std::string output() const
    {
      return getRoot()->output();
    }

    /*
     * Return true if all binary operators are '+' or '-'
     *
     */
    
    bool isAllSum();
    
  };
    
}

#endif
