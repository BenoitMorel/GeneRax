//
// File IoFormat.h
// Author : Julien Dutheil
// Created on: 2005
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#ifndef _IOFORMAT_H_
#define _IOFORMAT_H_

#include "../Exceptions.h"

// From STL:
#include <string>

namespace bpp
{

  /**
   * @brief The IOFormat interface.
   *
   * This is the most basal class of any format I/O implementation.
   */
  class IOFormat
  {
  public:
    IOFormat() {}
    virtual ~IOFormat() {}

  public:

    /**
     * @brief Get the type of data this format deals with.
     *
     * @return The type of data.
     */
    virtual const std::string getDataType() const = 0;
		
    /**
     * @brief Get the name of the file format.
     *
     * @return The name of the format implemented.
     */
    virtual const std::string getFormatName() const = 0;
		
    /**
     * @brief Get a description of the file format.
     *
     * @return A description of the format implemented.
     */
    virtual const std::string getFormatDescription() const = 0;
  };

} //end of namespace bpp.

#endif	// _IOFORMAT_H_

