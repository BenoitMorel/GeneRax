//
// File: ComputationTree.cpp
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

#include "ComputationTree.h"

#include "BinaryOperator.h"
#include "FunctionOperator.h"
#include "ConstantOperator.h"
#include "MathOperator.h"

#include <algorithm>

using namespace std;
using namespace bpp;

ComputationTree::ComputationTree(const std::string& formula, const std::map<std::string, Function*>& functionNames):
  AssociationTreeGlobalGraphObserver<Operator,short>(true)
{
  getGraph() ;

  std::string str2 = formula;
  
  str2.erase(std::remove_if(str2.begin(), 
                            str2.end(),
                            [](char x){return std::isspace(x);}),
             str2.end());

  setRoot(readFormula_(str2, functionNames));
}

std::shared_ptr<Operator> ComputationTree::readFormula_(const std::string& formula, const std::map<std::string, Function*>& functionNames)
{
  unsigned int level = 0;
  //inside parentheses check
  //case + or -
  //most right '+' or '-' (but not inside '()') search and split

  for(size_t i=formula.size();i>0;--i){

    char c = formula[i-1]; 

    if(c == ')'){
      ++level;
      continue;
    }
    
    if(c == '('){
      --level;
      continue;
    }

    if(level>0)
      continue;

    if ((c == '+' || c == '-') && !(i==1 || formula[i-2] == '*' || formula[i-2] == '/'
                                  || formula[i-2] == '+' || formula[i-2] == '-'))
    {

      //if sign

      std::shared_ptr<Operator> left=readFormula_(formula.substr(0,i-1), functionNames);
      std::shared_ptr<Operator> right=readFormula_(formula.substr(i), functionNames);

      shared_ptr<Operator> here(new BinaryOperator(c,left,right));
      
      createNode(here);
      
      setFather(left, here);
      setFather(right, here);
      
      return here;
    }
  }

  //case * or /
  //most right '*' or '/' (but not inside '()') search and split
  for(size_t i=formula.size();i>0;--i){
     
    char c = formula[i-1];

    if(c == ')'){
      ++level;
      continue;
    }
    
    if(c == '('){
      --level;
      continue;
    }

    if(level>0)
      continue;

    if(c == '*' || c == '/'){

      std::shared_ptr<Operator> left=readFormula_(formula.substr(0,i-1), functionNames);
      std::shared_ptr<Operator> right=readFormula_(formula.substr(i), functionNames);

      shared_ptr<Operator> here(new BinaryOperator(c,left,right));
      
      createNode(here);
      
      setFather(left, here);
      setFather(right, here);
      
      return here;
    }

  }


  if(formula[0]=='(')
    return readFormula_(formula.substr(1, formula.size()-2), functionNames);
  
  else
    //case value
  {
    shared_ptr<Operator> here;
    try
    {
      double v = TextTools::toDouble(formula);
      here = shared_ptr<Operator>(new ConstantOperator(v));
    }
    catch (Exception e)
    {
      std::map<std::string, Function*>::const_iterator it(functionNames.find(formula));

      if (it!=functionNames.end())
      {
        if (dynamic_cast<const DerivableSecondOrder*>(it->second))
          here=shared_ptr<Operator>(new FunctionOperator<DerivableSecondOrder>(*dynamic_cast<DerivableSecondOrder*>(it->second),formula));
        else
          if (dynamic_cast<const DerivableFirstOrder*>(it->second))
            here=shared_ptr<Operator>(new FunctionOperator<DerivableFirstOrder>(*dynamic_cast<DerivableFirstOrder*>(it->second),formula));
        else
          here=shared_ptr<Operator>(new FunctionOperator<Function>(*it->second,formula));
      }
      else
      {
        size_t posp=formula.find("(");
        if (posp==string::npos)
          throw Exception("ComputationTree::readFormula_ : unknown formula : "+ formula);
        
        std::shared_ptr<Operator> son=readFormula_(formula.substr(posp), functionNames);
        string fonc=formula.substr(0,posp);
        
        if (fonc=="exp")
          here=shared_ptr<Operator>(new MathOperator(&exp,"exp",son));
        else if (fonc=="log")
          here=shared_ptr<Operator>(new MathOperator(&log,"log",son));
        else
          throw Exception("ComputationTree::readFormula_ : unknown formula : "+ formula);
      }
    }
    
    this->getGraph();
    createNode(here);
    return here;
  }
  
  return NULL;
  //never
}

bool ComputationTree::isAllSum()
{
  std::unique_ptr<NodeIterator> it=allNodesIterator();

  for (;it->end();it->next())
  {
    const BinaryOperator* op = dynamic_cast<const BinaryOperator*>((**it).get());
    if (op && op->getSymbol()!='+' && op->getSymbol()!='-')
      return false;
  }
  return true;
  
}
