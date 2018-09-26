//
// File: BppVector.h
// Created by: Julien Dutheil
// Created on: Mon Apr 07 15:14 2008
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide utilitary
classes. This file belongs to the Bio++ Project.

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

#ifndef _BPPVECTOR_H_
#define _BPPVECTOR_H_

#include "Clonable.h"

// From the STL:
#include <vector>

namespace bpp
{

/**
 * @brief The BppVector object class.
 *
 * This class extends the std::vector class to support the Clonable interface.
 */
template<class TYPE>
class BppVector: public Clonable
{
  private:
    std::vector<TYPE> vector_;

  public:

    BppVector():
      vector_() {}
    
    BppVector(typename std::vector<TYPE>::size_type num, const TYPE& val = TYPE() ):
      vector_(num, val) {}

    BppVector(typename std::vector<TYPE>::iterator start, typename std::vector<TYPE>::iterator end):
      vector_(start, end) {}

    virtual ~BppVector() {}
  
  public:
  
    /**
     * @name The Clonable interface.
     *
     * @{
     */
    BppVector<TYPE> * clone() const { return new BppVector<TYPE>(*this); }
    /** @} */

    const std::vector<TYPE>& toSTL() const { return vector_; }
    
    std::vector<TYPE>& toSTL() { return vector_; }

    unsigned int size() const { return vector_.size(); }

    TYPE& operator[] (unsigned int i) { return vector_[i]; }
    
    const TYPE& operator[] (unsigned int i) const { return vector_[i]; }
  
};

} //end of namespace bpp.

#endif  //_BPPVECTOR_H_

