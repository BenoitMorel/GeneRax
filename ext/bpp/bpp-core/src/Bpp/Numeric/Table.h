//
// File: Table.h
// Created by: Laurent Guéguen
// Created on: dimanche 2 avril 2017, à 22h 59
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

#ifndef _TABLE_H_
#define _TABLE_H_

#include "VectorTools.h"
#include "TableExceptions.h"
#include "../Clonable.h"

// From the STL:
#include <string>
#include <vector>
#include <map>
#include <memory>


namespace bpp
{
/**
 * @brief This class corresponds to a 'dataset', <i>i.e.</i> a table
 * with data by rows and variable by columns.
 *
 * Data are stored as T objects, by column.
 */

  template<class T>
  class Table :
    public Clonable
  {
  protected:
    size_t nRow_, nCol_;
    std::vector< std::vector<T> > data_;
    std::vector<std::string> rowNames_;
    std::vector<std::string> colNames_;

  public:
    /**
     * @brief Build a new void Table object with nRow rows and nCol columns.
     *
     * @param nRow The number of rows of the Table.
     * @param nCol The number of columns of the Table.
     */
    
    Table(size_t nRow, size_t nCol) :
      nRow_(nRow),
      nCol_(nCol),
      data_(nCol),
      rowNames_(),
      colNames_()
    {
      for (size_t i = 0; i < nCol; i++)
      {
        data_[i].resize(nRow);
      }
    }

    /**
     * @brief Build a new void Table object with named columns.
     *
     * @param colNames The names of the columns of the Table.
     * @throw DuplicatedTableColumnNameException If colnames contains identical names.
     */

    Table(const std::vector<std::string>& colNames) throw (DuplicatedTableColumnNameException) :
      nRow_(0),
      nCol_(colNames.size()),
      data_(colNames.size()),
      rowNames_(),
      colNames_()
    {
      setColumnNames(colNames); // May throw an exception.
    }
    

    Table(const Table& table) :
      nRow_(table.nRow_),
      nCol_(table.nCol_),
      data_(table.data_),
      rowNames_(table.rowNames_),
      colNames_(table.colNames_)
    {
    }

    Table(const std::vector<std::vector<T> >& vt) :
      nRow_(vt.size()==0?0:vt[0].size()),      
      nCol_(vt.size()),
      data_(vt),
      rowNames_(),
      colNames_()
    {}

    Table& operator=(const Table& table)
    {
      nRow_ = table.nRow_;
      nCol_ = table.nCol_;
      data_ = table.data_;
      rowNames_ = table.rowNames_;
      colNames_ = table.colNames_;
      return *this;
    }

    Table& operator=(const std::vector<std::vector<T> >& vt)
    {
      nCol_ = vt.size();
      if (vt.size()==0)
      {
        nRow_ = 0;
        data_.clear();
      }
      else
      {
        nRow_=vt[0].size();
        data_=vt;
      }
      
      rowNames_.clear();
      colNames_.clear();
      return *this;
    }

    Table* clone() const { return new Table(*this); }

    virtual ~Table() {}
    
  public:
    /**
     * @return The data.
     *
     */

    const std::vector< std::vector<T> >& getData() const
    {
      return data_;
    }
    

    /**
     * @return The element at a given position.
     * @param rowIndex Row number.
     * @param colIndex Column number.
     * @throw IndexOutOfBoundsException If one of the index is greater or equal to the corresponding number of columns/rows.
     */

    T& operator()(size_t rowIndex, size_t colIndex) throw (IndexOutOfBoundsException)
    {
      if (colIndex >= nCol_)
        throw IndexOutOfBoundsException("Table::operator(size_t, size_t).", colIndex, 0, nCol_ - 1);
      if (rowIndex >= data_[colIndex].size())
        throw IndexOutOfBoundsException("Table::operator(size_t, size_t).", rowIndex, 0, data_[colIndex].size() - 1);
      return data_[colIndex][rowIndex];
    }

    /**
     * @return The element at a given position.
     * @param rowIndex Row number.
     * @param colIndex Column number.
     * @throw IndexOutOfBoundsException If one of the index is greater or equal to the corresponding number of columns/rows.
     */

    const T& operator()(size_t rowIndex, size_t colIndex) const throw (IndexOutOfBoundsException)
    {
      if (colIndex >= nCol_)
        throw IndexOutOfBoundsException("Table::operator(size_t, size_t).", colIndex, 0, nCol_ - 1);
      if (rowIndex >= data_[colIndex].size())
        throw IndexOutOfBoundsException("Table::operator(size_t, size_t).", rowIndex, 0, data_[colIndex].size() - 1);
      return data_[colIndex][rowIndex];
    }

    /**
     * @return The element at a given position.
     * @param rowName Row name.
     * @param colName Column name.
     * @throw NoTableRowNamesException If the table does not have names associated to rows.
     * @throw NoTableColumnNamesException If the table does not have names associated to columns.
     * @throw TableNameNotFoundException If one of rowName or colName do not match existing names.
     */

    T& operator()(const std::string& rowName, const std::string& colName)
      throw (NoTableRowNamesException, NoTableColumnNamesException, TableNameNotFoundException)
    {
      if (rowNames_.size() == 0)
        throw NoTableRowNamesException("Table::operator(const string &, const string &).");
      if (colNames_.size() == 0)
        throw NoTableColumnNamesException("Table::operator(const string &, const string &).");
      try
      {
        size_t rowIndex = VectorTools::which(rowNames_, rowName);
        size_t colIndex = VectorTools::which(colNames_, colName);
        return (*this)(rowIndex, colIndex);
      }
      catch (ElementNotFoundException<std::string>& ex)
      {
        throw TableNameNotFoundException("Table::operator(const string &, const string &).", *ex.getElement());
      }
    }
    
    /**
     * @return The element at a given position.
     * @param rowName Row name.
     * @param colName Column name.
     * @throw NoTableRowNamesException If the table does not have names associated to rows.
     * @throw NoTableColumnNamesException If the table does not have names associated to columns.
     * @throw TableNameNotFoundException If one of rowName or colName do
     * not match existing names.
     */
    const T& operator()(const std::string& rowName, const std::string& colName) const
      throw (NoTableRowNamesException, NoTableColumnNamesException, TableNameNotFoundException)
    {
      if (rowNames_.size() == 0)
        throw NoTableRowNamesException("Table::operator(const string &, const string &).");
      if (colNames_.size() == 0)
        throw NoTableColumnNamesException("Table::operator(const string &, const string &).");
      try
      {
        size_t rowIndex = VectorTools::which(rowNames_, rowName);
        size_t colIndex = VectorTools::which(colNames_, colName);
        return (*this)(rowIndex, colIndex);
      }
      catch (ElementNotFoundException<std::string>& ex)
      {
        throw TableNameNotFoundException("Table::operator(const string &, const string &).", *ex.getElement());
      }
    }
  
    

    /**
     * @return The element at a given position.
     * @param rowName Row name.
     * @param colIndex Column number.
     * @throw NoTableRowNamesException If the table does not have names associated to rows.
     * @throw IndexOutOfBoundsException If the index is greater or equal to the number of columns.
     * @throw TableNameNotFoundException If rowName do not match existing names.
     */

    T& operator()(const std::string& rowName, size_t colIndex)
      throw (NoTableRowNamesException, TableNameNotFoundException, IndexOutOfBoundsException)
    {
      if (rowNames_.size() == 0)
        throw NoTableRowNamesException("Table::operator(const string &, size_t).");
      if (colIndex >= nCol_)
        throw IndexOutOfBoundsException("Table::operator(const string &, size_t).", colIndex, 0, nCol_ - 1);
      try
      {
        size_t rowIndex = VectorTools::which(rowNames_, rowName);
        return (*this)(rowIndex, colIndex);
      }
      catch (ElementNotFoundException<std::string>& ex)
      {
        throw TableNameNotFoundException("Table::operator(const string &, size_t).", *ex.getElement());
      }
    }


    /**
     * @return The element at a given position.
     * @param rowName Row name.
     * @param colIndex Column number.
     * @throw NoTableRowNamesException If the table does not have names associated to rows.
     * @throw IndexOutOfBoundsException If the index is greater or equal to the number of columns.
     * @throw TableNameNotFoundException If rowName do not match existing names.
     */

    const T& operator()(const std::string& rowName, size_t colIndex) const
      throw (NoTableRowNamesException, TableNameNotFoundException, IndexOutOfBoundsException)
    {
      if (rowNames_.size() == 0)
        throw NoTableRowNamesException("Table::operator(const string &, size_t).");
      if (colIndex >= nCol_)
        throw IndexOutOfBoundsException("Table::operator(const string &, size_t).", colIndex, 0, nCol_ - 1);
      try
      {
        size_t rowIndex = VectorTools::which(rowNames_, rowName);
        return (*this)(rowIndex, colIndex);
      }
      catch (ElementNotFoundException<std::string>& ex)
      {
        throw TableNameNotFoundException("Table::operator(const string &, size_t).", *ex.getElement());
      }
    }


    /**
     * @return The element at a given position.
     * @param rowIndex Row number.
     * @param colName Column name.
     * @throw IndexOutOfBoundsException If the index is greater or equal to the number of rows.
     * @throw NoTableColumnNamesException If the table does not have names associated to columns.
     * @throw TableNameNotFoundException If colName do not match existing names.
     */

    T& operator()(size_t rowIndex, const std::string& colName)
      throw (IndexOutOfBoundsException, NoTableColumnNamesException, TableNameNotFoundException)
    {
      if (colNames_.size() == 0)
        throw NoTableColumnNamesException("Table::operator(size_t, const string &).");
      try
      {
        size_t colIndex = VectorTools::which(colNames_, colName);
        if (rowIndex >= data_[colIndex].size())
          throw IndexOutOfBoundsException("Table::operator(size_t, const string &).", rowIndex, 0, data_[colIndex].size() - 1);
        return (*this)(rowIndex, colIndex);
      }
      catch (ElementNotFoundException<std::string>& ex)
      {
        throw TableNameNotFoundException("Table::operator(const string &, const string &).", *ex.getElement());
      }
    }


    /**
     * @return The element at a given position.
     * @param rowIndex Row number.
     * @param colName Column name.
     * @throw IndexOutOfBoundsException If the index is greater or equal to the number of rows.
     * @throw NoTableColumnNamesException If the table does not have names associated to columns.
     * @throw TableNameNotFoundException If colName do not match existing names.
     */

    const T& operator()(size_t rowIndex, const std::string& colName) const
      throw (IndexOutOfBoundsException, NoTableColumnNamesException, TableNameNotFoundException)
    {
      if (colNames_.size() == 0)
        throw NoTableColumnNamesException("Table::operator(size_t, const string &).");
      try
      {
        size_t colIndex = VectorTools::which(colNames_, colName);
        if (rowIndex >= data_[colIndex].size())
          throw IndexOutOfBoundsException("Table::operator(size_t, const string &).", rowIndex, 0, data_[colIndex].size() - 1);
        return (*this)(rowIndex, colIndex);
      }
      catch (ElementNotFoundException<std::string>& ex)
      {
        throw TableNameNotFoundException("Table::operator(const string &, const string &).", *ex.getElement());
      }
    }

    /**
     * @name Work on columns.
     *
     * @{
     */

    /**
     * @return The number of columns in this table.
     */
    size_t getNumberOfColumns() const { return nCol_; }

    /**
     * @brief Set the column names of this table.
     *
     * @param colNames The row names.
     * @throw DimensionException If the number of names do not match the number of columns in the table.
     * @throw DuplicatedTableColumnNameException If names are not unique.
     */
    void setColumnNames(const std::vector<std::string>& colNames) throw (DimensionException, DuplicatedTableColumnNameException)
    {
      if (!VectorTools::isUnique(colNames))
        throw DuplicatedTableColumnNameException("Table::setColumnNames(...). Column names must be unique.");
      if (colNames.size() != nCol_)
        throw DimensionException("Table::setColumnNames.", colNames.size(), nCol_);
      else
        
        colNames_ = colNames;
    }

    /**
     * @brief Get the column names of this table.
     *
     * @return The column names of this table.
     * @throw NoTableColumnNamesException If no column names are associated to this table.
     */

    const std::vector<std::string>& getColumnNames() const throw (NoTableColumnNamesException)
    {
      if (colNames_.size() == 0)
        throw NoTableColumnNamesException("Table::getColumnNames().");
      return colNames_;
    }

    std::vector<std::string>& getColumnNames() throw (NoTableColumnNamesException)
    {
      if (colNames_.size() == 0)
        throw NoTableColumnNamesException("Table::getColumnNames().");
      return colNames_;
    }

    /**
     * @brief Get a given column name.
     *
     * @param index The index of the column.
     * @return The column name associated to the given column.
     * @throw NoTableColumnNamesException If no column names are associated to this table.
     * @throw IndexOutOfBoundsException If index is >= number of columns.
     */
    std::string getColumnName(size_t index) const throw (NoTableColumnNamesException, IndexOutOfBoundsException)
    {
      if (colNames_.size() ==0)
        throw NoTableColumnNamesException("Table::getColumnName(size_t).");
      if (index >= nCol_)
        throw IndexOutOfBoundsException("Table::getColumnName(size_t).", index, 0, nCol_ - 1);
      return colNames_[index];
    }


    /**
     * @return true If column names are associated to this table.
     */
    bool hasColumnNames() const { return colNames_.size() != 0; }

    /**
     * @return The values in the given column.
     * @param index The index of the column.
     * @throw IndexOutOfBoundsException If index is >= number of columns.
     */
    std::vector<T>& getColumn(size_t index)
    {
      if (index >= nCol_)
        throw IndexOutOfBoundsException("Table::getColumn(size_t).", index, 0, nCol_ - 1);
      return data_[index];
    }

    /**
     * @return The values in the given column.
     * @param index The index of the column.
     * @throw IndexOutOfBoundsException If index is >= number of columns.
     */
    const std::vector<T>& getColumn(size_t index) const
    {
      if (index >= nCol_)
        throw IndexOutOfBoundsException("Table::getColumn(size_t).", index, 0, nCol_ - 1);
      return data_[index];
    }


    /**
     * @return The values in the given column.
     * @param colName The name of the column.
     * @throw NoTableColumnNamesException If no column names are associated to this table.
     * @throw TableColumnNameNotFoundException If colName do not match existing column names.
     */
    std::vector<T>& getColumn(const std::string& colName) throw (NoTableColumnNamesException, TableColumnNameNotFoundException)
    {
      if (colNames_.size() == 0)
        throw NoTableColumnNamesException("Table::getColumn(const string &).");
      try
      {
        size_t colIndex = VectorTools::which(colNames_, colName);
        return data_[colIndex];
      }
      catch (ElementNotFoundException<std::string>& ex)
      {
        throw TableColumnNameNotFoundException("Table::getColumn(const string &).", colName);
      }
    }

    /**
     * @return The values in the given column.
     * @param colName The name of the column.
     * @throw NoTableColumnNamesException If no column names are associated to this table.
     * @throw TableColumnNameNotFoundException If colName do not match existing column names.
     */
    const std::vector<T>& getColumn(const std::string& colName) const throw (NoTableColumnNamesException, TableColumnNameNotFoundException)
    {
      if (colNames_.size() == 0)
        throw NoTableColumnNamesException("Table::getColumn(const string &).");
      try
      {
        size_t colIndex = VectorTools::which(colNames_, colName);
        return data_[colIndex];
      }
      catch (ElementNotFoundException<std::string>& ex)
      {
        throw TableColumnNameNotFoundException("Table::getColumn(const string &).", colName);
      }
    }


    /**
     * @brief Tell is a given column exists.
     *
     * @param colName The name of the column to look for.
     * @return true if the column was found, false if not or if there are no column names.
     */
    bool hasColumn(const std::string& colName) const
    {
      if (colNames_.size() == 0)
        return false;
      for (size_t i = 0; i < colNames_.size(); i++)
      {
        if ((colNames_)[i] == colName)
          return true;
      }
      return false;
    }


    /**
     * @brief Delete the given column.
     *
     * @param index The index of the column.
     * @throw IndexOutOfBoundsException If index is >= number of columns.
     */
    void deleteColumn(size_t index) throw (IndexOutOfBoundsException)
    {
      if (index >= nCol_)
        throw IndexOutOfBoundsException("Table::deleteColumn(size_t).", index, 0, nCol_ - 1);
      data_.erase(data_.begin() + (size_t)index); 
      if (colNames_.size()!=0)
        colNames_.erase(colNames_.begin() + (size_t)(index));
      nCol_--;
    }

/**
 * @brief Delete the given column.
 *
 * @param colName The name of the column.
 * @throw NoTableColumnNamesException If no column names are associated to this table.
 * @throw TableColumnNameNotFoundException If colName do not match existing column names.
 */
    void deleteColumn(const std::string& colName) throw (NoTableColumnNamesException, TableColumnNameNotFoundException)
    {
      if ((colNames_.size()==0))
        throw NoTableColumnNamesException("Table::deleteColumn(const string &).");
      try
      {
        size_t colIndex = VectorTools::which(colNames_, colName);
        data_.erase(data_.begin() + (size_t)(colIndex));
        colNames_.erase(colNames_.begin() + (size_t)(colIndex));
        nCol_--;
      }
      catch (ElementNotFoundException<std::string>& ex)
      {
        throw TableColumnNameNotFoundException("Table::deleteColumn(const string &).", colName);
      }
    }


    /**
     * @brief Add a new column.
     *
     * @param newColumn The new column values.
     * @param pos       The position optional (default : -1 for nb columns)
     * @throw DimensionException If the number of values does not match the number of rows.
     * @throw TableColumnNamesException If the table has row names.
     */

    void addColumn(const std::vector<T>& newColumn, int pos=-1) throw (DimensionException, TableColumnNamesException)
    {
      if (pos>(int)nCol_)
        throw DimensionException("Table::addColumn.", pos, nCol_);
      if (pos==-1)
        pos=nCol_;

      if (colNames_.size())
        throw TableColumnNamesException("Table::addColumn. Table has column names.");
      if (newColumn.size() != nRow_)
        throw DimensionException("Table::addColumn.", newColumn.size(), nRow_);

      data_.insert(data_.begin()+pos, newColumn);
      nCol_++;
    }

    /**
     * @brief Add a new column.
     *
     * @param colName   The name of the column.
     * @param newColumn The new column values.
     * @param pos       The position optional (default : -1 for nb columns)
     * @throw DimensionException If the number of values does not match the number of rows.
     * @throw NoTableColumnNamesException If the table does not have row names.
     * @throw DuplicatedTableColumnNameException If colName is already used.
     */

    void addColumn(const std::string& colName, const std::vector<T>& newColumn, int pos=-1) throw (DimensionException, NoTableColumnNamesException, DuplicatedTableColumnNameException)
    {
      if (pos>(int)nCol_)
        throw DimensionException("Table::addColumn.", pos, nCol_);
      if (pos==-1)
        pos=nCol_;
      
      if ((colNames_.size()==0))
      {
        if (nCol_ == 0)
          colNames_ = new std::vector<std::string>();
        else
          throw NoTableColumnNamesException("Table::addColumn. Table has column names.");
      }
      if (newColumn.size() != nRow_)
        throw DimensionException("Table::addColumn.", newColumn.size(), nRow_);
      if (nCol_ > 0 && find(colNames_.begin(), colNames_.end(), colName) != colNames_.end())
        throw DuplicatedTableColumnNameException("Table::addColumn(const string &, const std::vector<string> &). Column names must be unique.");

      colNames_.insert(colNames_.begin()+pos,colName);
      data_.insert(data_.begin()+pos, newColumn);
      nCol_++;
    }

    /**
     * @brief Add a new column
     *
     * @param stream the input stream
     * @param sep The row delimiter
     * @param pos       The position optional (default : nb cols)
     * @param rowCol the indice of row where colnames are, starting
     * from 0 (default -1 means no such column)
     */
    
    void addColumn(std::string& st, const std::string& sep = "\t", int pos=-1, int rowCol = -1)
    {
      if (pos>(int)nCol_)
        throw DimensionException("Table::addColumn.", pos, nCol_);
      if (pos==-1)
        pos=nCol_;

      StringTokenizer stok(st, sep, false, true);
      std::vector<std::string> row(stok.getTokens().begin(), stok.getTokens().end());

      if (row.size()!=nRow_+(rowCol>=0)?1:0)
        throw BadIntegerException("Table::addColumn. Bad number of rows: ", row.size());

      std::vector<T> newColumn;
      
      for (size_t i=0; i<row.size(); i++)
      {
        if (i==rowCol)
        {
          std::string colName=row[i];
          if (find(colNames_.begin(), colNames_.end(), colName) != colNames_.end())
            throw DuplicatedTableColumnNameException("Table::addColumn(const std::vector<string> &). Column names must be unique.");

          colNames_.insert(colNames_.begin()+(size_t)pos,colName);
        }
        else
        {
          std::stringstream ss(row[i]);
          T t;
          ss >> t;
          newColumn.push_back(t);
        }
      }
      data_.insert(data_.begin()+pos, newColumn);

      nCol_++;      
    }


    /**
     * @brief Set a new column.
     *
     * @param newColumn The new column values.
     * @param pos       The position
     * @throw DimensionException If the number of values does not match the number of rows.
     * @throw TableColumnNamesException If the table has row names.
     */
    
    void setColumn(const std::vector<T>& newColumn, size_t pos) throw (DimensionException, TableColumnNamesException)
    {
      if (pos>=(int)nCol_)
        throw DimensionException("Table::setColumn.", pos, nCol_);

      if (newColumn.size() != nRow_)
        throw DimensionException("Table::setColumn.", newColumn.size(), nRow_);

      data_[pos]=newColumn;
    }

    /**
     * @brief Set a new column.
     *
     * @param colName   The name of the column.
     * @param newColumn The new column values.
     * @param pos       The position
     * @throw DimensionException If the number of values does not match the number of rows.
     * @throw NoTableColumnNamesException If the table does not have row names.
     * @throw DuplicatedTableColumnNameException If colName is already used.
     */
    void setColumn(const std::string& colName, const std::vector<T>& newColumn, size_t pos) throw (DimensionException, NoTableColumnNamesException, DuplicatedTableColumnNameException)
    {
      if (pos>=nCol_)
        throw DimensionException("Table::setColumn.", pos, nCol_);
      
      if (colNames_.size()==0)
        throw NoTableColumnNamesException("Table::setColumn. Table has column names.");

      if (newColumn.size() != nRow_)
        throw DimensionException("Table::setColumn.", newColumn.size(), nRow_);
      
      if (find(colNames_.begin(), colNames_.end(), colName) != colNames_.end() && 
          find(colNames_.begin(), colNames_.end(), colName) != colNames_.begin()+pos)
        throw DuplicatedTableColumnNameException("Table::setColumn(const string &, const std::vector<T> &, size_t). Column names must be unique.");

      colNames_[pos]=colName;
      data_[pos]=newColumn;
    }

    /** @} */

    /**
     * @name Work on rows.
     *
     * @{
     */

    /**
     * @return The number of rows in this table.
     */
    size_t getNumberOfRows() const { return nRow_; }

    /**
     * @brief Set the row names of this table.
     *
     * @param rowNames The row names.
     * @throw DimensionException If the number of names do not match the number of rows in the table.
     * @throw DuplicatedTableRowNameException If names are not unique.
     */
    void setRowNames(const std::vector<std::string>& rowNames) throw (DimensionException, DuplicatedTableRowNameException)
    { 
      if (!VectorTools::isUnique(rowNames))
      {
        throw DuplicatedTableRowNameException("Table::setRowNames(...). Row names must be unique.");
      }
      if (rowNames.size() != nRow_)
        throw DimensionException("Table::setRowNames.", rowNames.size(), nRow_);
      else
      {
        rowNames_ = rowNames;
      }
    }


    /**
     * @brief Get the row names of this table.
     *
     * @return The row names of this table.
     * @throw NoTableRowNamesException If no row names are associated to this table.
     */
    
    std::vector<std::string> getRowNames() throw (NoTableRowNamesException)
    {
      if (rowNames_.size() == 0)
        throw NoTableRowNamesException("Table::getRowNames().");
      return rowNames_;
    }

    const std::vector<std::string>& getRowNames() const throw (NoTableRowNamesException)
    {
      if (rowNames_.size() == 0)
        throw NoTableRowNamesException("Table::getRowNames().");
      return rowNames_;
    }


    /**
     * @brief Tell is a given row exists.
     *
     * @param rowName The name of the row to look for.
     * @return true if the row was found, false if not or if there are no row names.
     */
    bool hasRow(const std::string& rowName) const
    {
      if (rowNames_.size() == 0)
        return false;
      for (size_t i = 0; i < rowNames_.size(); i++)
      {
        if ((rowNames_)[i] == rowName)
          return true;
      }
      return false;
    }


    /**
     * @brief Get a given row name.
     *
     * @param index The index of the row.
     * @return The row name associated to the given row.
     * @throw NoTableRowNamesException If no row names are associated to this table.
     * @throw IndexOutOfBoundsException If index is >= number of rows.
     */
    std::string getRowName(size_t index) const throw (NoTableRowNamesException, IndexOutOfBoundsException)
    {
      if (rowNames_.size() == 0)
        throw NoTableRowNamesException("Table::getRowName(size_t).");
      if (index >= nRow_)
        throw IndexOutOfBoundsException("Table::getRowName(size_t).", index, 0, nRow_ - 1);
      return (rowNames_)[index];
    }


    /**
     * @return true If row names are associated to this table.
     */
    bool hasRowNames() const { return rowNames_.size() != 0; }

    /**
     * @return A vector which contains a copy  in the given row.
     * @param index The index of the row.
     * @throw IndexOutOfBoundsException If index is >= number of rows.
     */
    std::vector<T> getRow(size_t index) const throw (IndexOutOfBoundsException)
    {
      if (index >= nRow_)
        throw IndexOutOfBoundsException("Table::getRow(size_t).", index, 0, nRow_ - 1);
      std::vector<T> row;
      for (size_t i = 0; i < nCol_; i++)
      {
        row.push_back(data_[i][index]);
      }
      return row;
    }


    /**
     * @return A vector which contains a copy  in the given row.
     * @param rowName The name of the row.
     * @throw NoTableRowNamesException If no row names are associated to this table.
     * @throw TableRowNameNotFoundException If rowName do not match existing row names.
     */
    std::vector<T> getRow(const std::string& rowName) const throw (NoTableRowNamesException, TableRowNameNotFoundException)
    {
      if ((rowNames_.size()==0))
        throw NoTableRowNamesException("Table::getRow(const string &).");
      try
      {
        size_t rowIndex = VectorTools::which(rowNames_, rowName);
        std::vector<T> row;
        for (size_t i = 0; i < nCol_; i++)
        {
          row.push_back(data_[i][rowIndex]);
        }
        return row;
      }
      catch (ElementNotFoundException<std::string>& ex)
      {
        throw TableRowNameNotFoundException("Table::getRow(const string &).", rowName);
      }
    }


    /**
     * @brief Delete the given row.
     *
     * @param index The index of the row.
     * @throw IndexOutOfBoundsException If index is >= number of row.
     */
    void deleteRow(size_t index) throw (IndexOutOfBoundsException)
    {
      for (size_t j = 0; j < nCol_; j++)
      {
        std::vector<T>* column = &data_[j];
        if (index >= column->size())
          throw IndexOutOfBoundsException("Table::deleteRow(size_t).", index, 0, column->size() - 1);
        column->erase(column->begin() + index);
      }
      if (rowNames_.size()!=0)
        rowNames_.erase(rowNames_.begin() + index);
      nRow_--;
    }

    /**
     * @brief Delete the given rows.
     *
     * @param index The index of the first row.
     * @param len the number of rows to delete
     * @throw IndexOutOfBoundsException If index is >= number of row.
     */
    
    void deleteRows(size_t index, size_t len) throw (IndexOutOfBoundsException)
    {
      for (size_t j = 0; j < nCol_; j++)
      {
        std::vector<T>* column = &data_[j];
        if (index >= column->size())
          throw IndexOutOfBoundsException("Table::deleteRow(size_t).", index, 0, column->size() - 1);
        column->erase(column->begin() + index, column->begin() + index + len);
      }
      if (rowNames_.size()!=0)
        rowNames_.erase(rowNames_.begin() + index, rowNames_.begin() + index + len);
      nRow_--;
    }


    /**
     * @brief Delete the given row.
     *
     * @param rowName The name of the row.
     * @throw NoTableRowNamesException If no row names are associated to this table.
     * @throw TableRowNameNotFoundException If rowName do not match existing column names.
     */
    void deleteRow(const std::string& rowName) throw (NoTableRowNamesException, TableRowNameNotFoundException)
    {
      if ((rowNames_.size()==0))
        throw NoTableRowNamesException("Table::deleteRow(const string &).");
      try
      {
        size_t rowIndex = VectorTools::which(rowNames_, rowName);
        for (size_t j = 0; j < nCol_; j++)
        {
          std::vector<T>* column = &data_[j];
          column->erase(column->begin() + rowIndex);
        }
        rowNames_.erase(rowNames_.begin() + rowIndex);
        nRow_--;
      }
      catch (ElementNotFoundException<std::string>& ex)
      {
        throw TableRowNameNotFoundException("Table::deleteRow(const string &).", rowName);
      }
    }


    /**
     * @brief Add a new row.
     *
     * @param newRow The new row values.
     * @param pos       The position optional (default : -1 for nb rows)
     */
    void addRow(const std::vector<T>& newRow, int pos=-1)
    {
      if (pos>(int)nRow_)
        throw DimensionException("Table::addRow.", pos, nRow_);

      if (pos==-1)
        pos=nRow_;
      
      if (rowNames_.size()!=0)
        throw TableRowNamesException("Table::addRow. Table has row names.");
      if (newRow.size() != nCol_)
        throw DimensionException("Table::addRow.", newRow.size(), nCol_);
      for (size_t j = 0; j < nCol_; j++)
      {
        data_[j].insert(data_[j].begin()+(size_t)pos, newRow[j]);
      }
      nRow_++;
    }

    /**
     * @brief Add a new row at a given position
     *
     * @param rowName   The name of the row.
     * @param newRow    The new row values.
     * @param pos       The position optional (default : -1 for nb rows)
     */

    void addRow(const std::string& rowName, const std::vector<T>& newRow, int pos=-1)
    {
      if (pos>(int)nRow_)
        throw DimensionException("Table::addRow.", pos, nRow_);
      if (pos==-1)
        pos=nRow_;

      if ((rowNames_.size()==0))
      {
        if (nRow_ == 0)
          rowNames_ = new std::vector<std::string>();
        else
          throw NoTableRowNamesException("Table::addRow. Table has row names.");
      }
      if (newRow.size() != nCol_)
        throw DimensionException("Table::addRow.", newRow.size(), nCol_);
      if (nRow_ > 0 && find(rowNames_.begin(), rowNames_.end(), rowName) != rowNames_.end())
        throw DuplicatedTableRowNameException("Table::addRow(const string &, const std::vector<string> &). Row names must be unique.");
      rowNames_.insert(rowNames_.begin()+(size_t)pos,rowName);
      for (size_t j = 0; j < nCol_; j++)
      {
        data_[j].insert(data_[j].begin()+(size_t)pos, newRow[j]);
      }
      nRow_++;
    }

    /**
     * @brief Add a new row.
     *
     * @param stream the input stream
     * @param sep The column delimiter
     * @param pos       The position optional (default : nb rows)
     * @param rowCol the indice of column where rownames are, starting
     * from 0 (default -1 means no such column)
     */
    
    void addRow(std::string& st, const std::string& sep = "\t", int pos=-1, int rowCol = -1)
    {
      if (pos>(int)nRow_)
        throw DimensionException("Table::addRow.", pos, nRow_);
      if (pos==-1)
        pos=nRow_;

      StringTokenizer stok(st, sep, false, true);
      std::vector<std::string> row(stok.getTokens().begin(), stok.getTokens().end());

      if (row.size()!=nCol_+(rowCol>=0)?1:0)
        throw BadIntegerException("Table::addRow. Bad number of columns: ", row.size());

      size_t id=0;
      
      for (size_t i=0; i<row.size(); i++)
      {
        if (i==rowCol)
        {
          std::string rowName=row[i];
          if (find(rowNames_.begin(), rowNames_.end(), rowName) != rowNames_.end())
            throw DuplicatedTableRowNameException("Table::addRow(const std::vector<string> &). Row names must be unique.");

          rowNames_.insert(rowNames_.begin()+(size_t)pos,rowName);
        }
        else
        {
          std::stringstream ss(row[i]);
          T t;
          ss >> t;
          data_[id].insert(data_[id].begin()+(size_t)pos,t);
          id++;
        }
      }
      nRow_++;      
    }


    /** @} */

  public:
    /**
     * @brief Read a table form a stream in CSV-like
     *
     * By default, if the first line as one column less than the second one,
     * the first line is taken as column names, and the first column as row names.
     * Otherwise, no column names and no row names are specified, unless
     * explicitely precised by the user.
     * Note: if byRow=false, exchange row and column roles.
     *
     * @param in       The input stream.
     * @param byRow    Tell if the table is filled by rows
     * @param sep      The line delimiter.
     * @param header   Tell if the first line must be used for names,
     * otherwise use default.
     * @param names     Use a column as rowNames (or a row as columnNames). If positive, use the specified column to compute rownames, otherwise use default;
     * @return         A pointer toward a new Table object.
     */
    
    static Table<T>* read(std::istream& in, bool byRow, const std::string& sep = "\t", bool header = true, int names = -1)
    {
      std::string firstLine  = FileTools::getNextLine(in);
      StringTokenizer st1(firstLine, sep, false, true);
      std::vector<std::string> row1(st1.getTokens().begin(), st1.getTokens().end());
      size_t nCol = row1.size();
      Table<T>* dt;

      if (header)
      {
        // Use first line as header.
        if (byRow)
        {
          dt = new Table<T>(0,nCol);
          dt->setColumnNames(row1);
        }
        else
        {
          dt = new Table<T>(nCol,0);
          dt->setRowNames(row1);
        }
      }
      else
      {
        if (byRow)
        {
          dt = new Table<T>(0,nCol-(names>=0?1:0));
          dt->addRow(firstLine, sep, names);
        }
        else
        {
          dt = new Table<T>(nCol-(names>=0?1:0),0);
          dt->addColumn(firstLine, sep, names);
        }
      }
            
      // Now read each line:
      std::string line = FileTools::getNextLine(in);
      
      while (!TextTools::isEmpty(line))
      {
        if (byRow)
          dt->addRow(line, sep, names);
        else
          dt->addColumn(line, sep, names);
          
        line = FileTools::getNextLine(in);
      }
      return dt;
    }


    /**
     * @brief Write a Table object to stream in CVS-like format.
     *
     * @param data         The table to write.
     * @param out          The output stream.
     * @param byRow        Tell if the table is written by rows.
     * @param sep          The column delimiter.
     * @param alignHeaders If true, add a delimiter before the first column header if there is row names.
     */
    static void write(const Table& data, std::ostream& out, bool byRow, const std::string& sep = "\t", bool alignHeaders = false)
    {
      size_t n = (byRow?data.getNumberOfColumns():data.getNumberOfRows());
      if (n == 0)
        return;
      size_t m = (byRow?data.getNumberOfRows():data.getNumberOfColumns());

      bool frontNames=((byRow && data.hasRowNames()) || (!byRow && data.hasColumnNames()));

      if (byRow && data.hasColumnNames() || (!byRow && data.hasRowNames()))
      { // Write header
        std::vector<std::string> names = (byRow?data.getColumnNames():data.getRowNames());
        
        if (alignHeaders && frontNames)
          out << sep;
        out << names[0];
        for (size_t i = 1; i < n; i++)
        {
          out << sep << names[i];
        }
        out << std::endl;
      }
      
      // now write each row (or column):
        
      for (size_t i = 0; i < m; i++)
      {
        if (frontNames){
          out << (byRow?data.getRowName(i):data.getColumnName(i)) << sep;
        }

        out << (byRow?data(i, 0):data(0, i));
        for (size_t j = 1; j < n; j++)
        {
          out << sep << (byRow?data(i, j):data(j, i));
        }
        out << std::endl;
      }
    }

    static void write(const Table& data, bpp::OutputStream& out, bool byRow, const std::string& sep = "\t", bool alignHeaders = false)
    {
      size_t n = (byRow?data.getNumberOfColumns():data.getNumberOfRows());
      if (n == 0)
        return;
      size_t m = (byRow?data.getNumberOfRows():data.getNumberOfColumns());

      bool frontNames=((byRow && data.hasRowNames()) || (!byRow && data.hasColumnNames()));

      if (byRow && data.hasColumnNames() || (!byRow && data.hasRowNames()))
      { // Write header
        std::vector<std::string> names = (byRow?data.getColumnNames():data.getRowNames());
        
        if (alignHeaders && frontNames)
          out << sep;
        out << names[0];
        for (size_t i = 1; i < n; i++)
        {
          out << sep << names[i];
        }
        out.endLine();
      }
      
      // now write each row (or column):

      for (size_t i = 0; i < m; i++)
      {
        if (frontNames)
          out << (byRow?data.getRowName(i):data.getColumnName(i) << sep);
        
        out << (byRow?data(i, 0):data(0, i));
        for (size_t j = 1; j < n; j++)
        {
          out << sep << byRow?data(i, j):data(j, i);
        }
        out.endLine();
      }
    }

  };
} // end of namespace bpp.

#endif // _Table_H_

