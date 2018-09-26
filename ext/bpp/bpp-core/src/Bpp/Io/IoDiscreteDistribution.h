//
// File: IoDiscreteDistribution.h
// Created by: Laurent Guéguen
// Created on: lundi 3 septembre 2012, à 14h 35
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

#ifndef _IODISCRETEDISTRIBUTION_H_
#define _IODISCRETEDISTRIBUTION_H_

#include "../Numeric/Prob/DiscreteDistribution.h"
#include "../Exceptions.h"
#include "IoFormat.h"
#include "OutputStream.h"

namespace bpp
{
  /**
   * @brief General interface for model I/O.
   */
  class IoDiscreteDistribution:
    public virtual IOFormat
  {
  public:
    IoDiscreteDistribution() {}
    virtual ~IoDiscreteDistribution() {}

  public:
    virtual const std::string getDataType() const { return "Discrete Distribution"; }
  };

  /**
   * @brief General interface for distance matrix readers.
   */
  class IDiscreteDistribution:
    public virtual IoDiscreteDistribution
  {
  public:
    IDiscreteDistribution() {}
    virtual ~IDiscreteDistribution() {}

  public:
    /**
     * @brief Read a discrete distribution from a string.
     *
     * @param distrDescription A string describing the distribution in the format.
     * @param parseArguments Attempt to parse function arguments. If not, only store them and use default values instead.
     * @return A new DiscreteDistribution object according to options specified.
     * @throw Exception if an error occured.
     */
    virtual DiscreteDistribution* read(
        const std::string& distrDescription,
        bool parseArguments = true) = 0;

    /**
     * @return The arguments and their unparsed values from the last call of the read function, if there are any.
     */
    virtual const std::map<std::string, std::string>& getUnparsedArguments() const = 0;

  };

  /**
   * @brief General interface writers.
   */
  class ODiscreteDistribution:
    public virtual IoDiscreteDistribution
  {
  public:
    ODiscreteDistribution() {}
    virtual ~ODiscreteDistribution() {}

  public:
    /**
     * @brief Write a discrete distribution to a stream.
     *
     * @param dist A discrete distribution object;
     * @param out The output stream;
     * @param globalAliases parameters linked to global alias. 
     * @param writtenNames is the vector of the written
     *        parameters so far [in, out];
     * @throw Exception if an error occured.
     */
    virtual void write(const DiscreteDistribution& dist,
                       OutputStream& out,
                       std::map<std::string, std::string>& globalAliases,
                       std::vector<std::string>& writtenNames) const = 0;
  };


} //end of namespace bpp.

#endif //_IODISCRETEDISTRIBUTION_H_

