//
// File: DataTableExceptions.h
// Created by: Julien Dutheil
// Created on: Tue Nov 2005 14:10
// from file DataTable.h
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#ifndef _TABLEEXCEPTIONS_H_
#define _TABLEEXCEPTIONS_H_

//#include "VectorTools.h"

// From Utils:
#include "../Exceptions.h"
#include "../Text/TextTools.h"

// From the STL:
#include <string>

namespace bpp
{

/**
 * @brief Exception thrown when a given name is not found is a DataTable object.
 */
  class TableNameNotFoundException:
    public Exception
  {
  protected:
    std::string _name;
		
  public:
    TableNameNotFoundException(const std::string & text, const std::string & name) :
      Exception("TableNameNotFoundException: " + name + ". " + text), _name(name) {}
    virtual ~TableNameNotFoundException() throw() {}

  public:
    std::string getName() const { return _name; }		
  };

/**
 * @brief Exception thrown when a given row name is not found is a DataTable object.
 */
  class TableRowNameNotFoundException:
    public TableNameNotFoundException
  {
  public:
    TableRowNameNotFoundException(const std::string & text, const std::string & name) :
      TableNameNotFoundException("TableRowNameNotFoundException: " + name + ". " + text, name) {}
    virtual ~TableRowNameNotFoundException() throw() {}
  };

/**
 * @brief Exception thrown when a given column name is not found is a DataTable object.
 */
  class TableColumnNameNotFoundException:
    public TableNameNotFoundException
  {
  public:
    TableColumnNameNotFoundException(const std::string & text, const std::string & name) :
      TableNameNotFoundException("TableColumnNameNotFoundException: " + name + ". " + text, name) {}
    virtual ~TableColumnNameNotFoundException() throw() {}
  };

/**
 * @brief Exception thrown when trying to retrieve a row by its name
 * and no row names have been specified.
 */
  class NoTableRowNamesException:
    public Exception
  {
  public:
    NoTableRowNamesException(const std::string & text) :
      Exception("NoTableRowNamesException: "+text) {}
    virtual ~NoTableRowNamesException() throw() {}
  };

/**
 * @brief Exception thrown when trying to retrieve a column by its name
 * and no column names have been specified.
 */
  class NoTableColumnNamesException:
    public Exception
  {
  public:
    NoTableColumnNamesException(const std::string & text) :
      Exception("NoTableColumnNamesException: "+text) {}
    virtual ~NoTableColumnNamesException() throw() {}
  };

/**
 * @brief General exception class dealing with row names.
 */
  class TableRowNamesException:
    public Exception
  {
  public:
    TableRowNamesException(const std::string & text) :
      Exception("TableRowNamesException: "+text) {}
    virtual ~TableRowNamesException() throw() {}
  };

/**
 * @brief General exception class dealing with column names.
 */
  class TableColumnNamesException:
    public Exception
  {
  public:
    TableColumnNamesException(const std::string & text) :
      Exception("TableColumnNamesException: "+text) {}
    virtual ~TableColumnNamesException() throw() {}
  };

/**
 * @brief Exception thrown when attempting to duplicate a row name.
 */
  class DuplicatedTableRowNameException:
    public Exception
  {
  public:
    DuplicatedTableRowNameException(const std::string & text) :
      Exception("DuplicatedTableRowNameException: "+text) {}
    virtual ~DuplicatedTableRowNameException() throw() {}
  };

/**
 * @brief Excpetion thrown when attempting to duplicate a column name.
 */
  class DuplicatedTableColumnNameException:
    public Exception
  {
  public:
    DuplicatedTableColumnNameException(const std::string & text) :
      Exception("DuplicatedTableColumnNameException: "+text) {}
    virtual ~DuplicatedTableColumnNameException() throw() {}
  };

} //end of namespace bpp.

#endif //_TABLEEXCEPTIONS_H_

