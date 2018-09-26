//
// File: Bpp0ParametrizableFormat.h
// Created by: Laurent Guéguen
// Created on: lundi 3 septembre 2012, à 15h 30
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

#ifndef _BPPOPARAMETRIZABLEFORMAT_H_
#define _BPPOPARAMETRIZABLEFORMAT_H_

#include "IoParametrizable.h"

namespace bpp
{

  /**
   * @brief Parametrizable output in BppO format.
   *
   * Writes a new parametrizable object according to BppO description
   * syntax (see the Bio++ Progam Suite manual for a detailed
   * description of this syntax).
   *
   */
  class BppOParametrizableFormat:
    public OParametrizable
  {
  public:
    BppOParametrizableFormat() {}
    virtual ~BppOParametrizableFormat() {}

  public:
    const std::string getFormatName() const { return "BppO"; }

    const std::string getFormatDescription() const { return "Bpp Options format."; }

    /**
     * @brief Write a Parametrizable to a stream.
     *
     * @param parametrizable A pointer to a Parametrizable object;
     * @param out The output stream;
     * @param writtenNames is the vector of the written
     *        parameters so far [in, out];
     * @param printComma boolean if a comma should be written at the
     *        beginning of the description.
     */
    
    void write(const Parametrizable* parametrizable,
               OutputStream& out,
               std::vector<std::string>& writtenNames,
               bool printComma = false) const;
    
    /**
     * @brief Write a ParameterAliasable to a stream.
     *
     * @param parametrizable A pointer to a Parametrizable object;
     * @param out The output stream;
     * @param globalAliases parameters linked to global alias; 
     * @param names the names of the parameters to be written;
     * @param writtenNames is the vector of the written
     *        parameters so far [in, out];
     * @param printLocalAliases boolean if local aliases should be written;
     * @param printComma boolean if a comma should be written at the
     *        beginning of the description.
     */
    
    void write(const ParameterAliasable* parametrizable,
               OutputStream& out,
               std::map<std::string, std::string>& globalAliases,
               const std::vector<std::string>& names,
               std::vector<std::string>& writtenNames,
               bool printLocalAliases = true,
               bool printComma = false) const;
  };

} //end of namespace bpp.

#endif //_BPPOPARAMETRIZABLEFORMAT_H_

