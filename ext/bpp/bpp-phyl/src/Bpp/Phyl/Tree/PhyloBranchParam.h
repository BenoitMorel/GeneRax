//
// File: PhyloBranchParam.h
// Created by: Laurent Guéguen
// Created on: mercredi 14 septembre 2016, à 23h 29
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

#ifndef _PHYLOBRANCH_PARAM_H_
#define _PHYLOBRANCH_PARAM_H_

#include <Bpp/Numeric/AbstractParametrizable.h>
#include <Bpp/Numeric/NumConstants.h>


namespace bpp
{

  class PhyloBranch;

  /**
   * Basic branch in which length is parameterized.
   *
   * Note : Branch Length must be the first parameter (index 0).
   *
   **/
  
  class PhyloBranchParam:
    public AbstractParametrizable
  {
  public:
    /**
     * @brief Constructors.
     *
     * @warning phyloTree_ does not know the edge exists.
     *
     */
    
    PhyloBranchParam(const std::string& prefix = ""):
      AbstractParametrizable(prefix)
    {
      addParameter_(new Parameter("BrLen", NumConstants::SMALL(), &Parameter::R_PLUS_STAR, 0));
    }

    PhyloBranchParam(const PhyloBranch& branch);

    /**
     * @brief Copy constructor.
     *
     * @param branch The branch to copy.
     */
    PhyloBranchParam(const PhyloBranchParam& branch) :
      AbstractParametrizable(branch)
    {
    }
    
    
    /**
     * @brief Assignation operator.
     *
     * @warning This operator copies all fields, excepted father and
     * son branch pointers. Without specific id, the PhyloTree does
     * not know the existence of the new edge.
     *
     * @param branch the branch to copy.
     * @return A reference toward this branch.
     */
    PhyloBranchParam& operator=(const PhyloBranchParam& branch)
    {
      AbstractParametrizable::operator=(*this);
      return *this;
    }
    
    
    PhyloBranchParam* clone() const { return new PhyloBranchParam(*this); }
    
    
    /**
     * @brief What is the branch length?
     * @return a double representing the branch length, 0 if length is
     * not defined.
     *
     */

    double getLength() const
    {
      return getParameter_(0).getValue();
    }
    
    
    /**
     * Define a new branch length
     * @param newLength a double repserenting the new length of the branch
     */
    
    void setLength(double newLength)
    {
      getParameter_(0).setValue(newLength);
    }
    

    friend class ParametrizablePhyloTree;
    
    /** @} */
    
    
  }; //end of class PhyloBranchParam


} //end of namespace bpp.

#endif  //_PHYLOBRANCH_PARAM_H_
