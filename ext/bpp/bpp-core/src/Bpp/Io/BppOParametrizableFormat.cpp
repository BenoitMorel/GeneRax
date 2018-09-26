//
// File: BppOParametrizableFormat.cpp
// Created by: Laurent Guéguen
// Created on: lundi 3 septembre 2012, à 15h 37
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "BppOParametrizableFormat.h"

using namespace bpp;

// From the STL:
#include <iomanip>
#include <algorithm>


using namespace std;

void BppOParametrizableFormat::write(const Parametrizable* parametrizable,
                                     OutputStream& out,
                                     std::vector<std::string>& writtenNames,
                                     bool printComma) const
{
  ParameterList pl = parametrizable->getParameters();
  int p = out.getPrecision();
  out.setPrecision(12);
  bool flag = printComma;
  for (size_t i = 0; i < pl.size(); ++i)
  {
    if (find(writtenNames.begin(), writtenNames.end(), pl[i].getName()) == writtenNames.end())
    {
      if (flag)
        out << ",";
      else
        flag = true;
      string pname = parametrizable->getParameterNameWithoutNamespace(pl[i].getName());
        
      (out << pname << "=").enableScientificNotation(false) << pl[i].getValue();
        
    }
  }
  out.setPrecision(p);
}

void BppOParametrizableFormat::write(const ParameterAliasable* parametrizable,
                                     OutputStream& out,
                                     std::map<std::string, std::string>& globalAliases,
                                     const std::vector<std::string>& names,
                                     std::vector<std::string>& writtenNames,
                                     bool printLocalAliases,
                                     bool printComma) const
{
  ParameterList pl = parametrizable->getIndependentParameters().subList(names);
  int p = out.getPrecision();
  out.setPrecision(12);
  bool flag = printComma;
  for (size_t i = 0; i < pl.size(); ++i)
  {
    if (find(writtenNames.begin(), writtenNames.end(), pl[i].getName()) == writtenNames.end())
    {
      if (flag)
        out << ",";
      else
        flag = true;
      string pname = parametrizable->getParameterNameWithoutNamespace(pl[i].getName());
      
      // Check for global aliases:
      if (globalAliases.find(pl[i].getName()) == globalAliases.end())
      {
        (out << pname << "=").enableScientificNotation(false) << pl[i].getValue();
      }
      else
        out << pname << "=" << globalAliases[pl[i].getName()];
      
      // Now check for local aliases:
      if (printLocalAliases)
      {
        vector<string> aliases = parametrizable->getAlias(pname);
        for (size_t j = 0; j < aliases.size(); ++j)
        {
          out << ", " << aliases[j] << "=" << pname;
        }
      }
      writtenNames.push_back(pl[i].getName());
    }
  }
  out.setPrecision(p);
}
