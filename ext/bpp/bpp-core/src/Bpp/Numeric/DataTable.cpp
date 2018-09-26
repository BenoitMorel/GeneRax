//
// File: DataTable.cpp
// Created by: Julien Dutheil
// Created on: Aug 2005
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

#include "DataTable.h"
#include "VectorTools.h"
#include "../Io/FileTools.h"
#include "../Text/TextTools.h"
#include "../Text/StringTokenizer.h"

using namespace bpp;
using namespace std;

/******************************************************************************/

DataTable::DataTable(size_t nRow, size_t nCol) :
  nRow_(nRow),
  nCol_(nCol),
  data_(nCol),
  rowNames_(0),
  colNames_(0)
{
  for (size_t i = 0; i < nCol; i++)
  {
    data_[i].resize(nRow);
  }
}

DataTable::DataTable(size_t nCol) :
  nRow_(0),
  nCol_(nCol),
  data_(nCol),
  rowNames_(0),
  colNames_(0)
{}

DataTable::DataTable(const std::vector<std::string>& colNames) noexcept(false) :
  nRow_(0),
  nCol_(colNames.size()),
  data_(colNames.size()),
  rowNames_(0),
  colNames_(0)

{
  setColumnNames(colNames); // May throw an exception.
}

DataTable::DataTable(const DataTable& table) :
  nRow_(table.nRow_),
  nCol_(table.nCol_),
  data_(table.data_),
  rowNames_(0),
  colNames_(0)
{
  if (table.rowNames_)
    rowNames_ = new vector<string>(*table.rowNames_);
  if (table.colNames_)
    colNames_ = new vector<string>(*table.colNames_);
}

DataTable& DataTable::operator=(const DataTable& table)
{
  nRow_ = table.nRow_;
  nCol_ = table.nCol_;
  data_ = table.data_;
  if (rowNames_)
    delete rowNames_;
  if (colNames_)
    delete colNames_;
  rowNames_ = 0;
  colNames_ = 0;
  if (table.rowNames_)
    rowNames_ = new vector<string>(*table.rowNames_);
  if (table.colNames_)
    colNames_ = new vector<string>(*table.colNames_);
  return *this;
}

/******************************************************************************/

DataTable::~DataTable()
{
  if (rowNames_ != NULL)
    delete rowNames_;
  if (colNames_ != NULL)
    delete colNames_;
}

/******************************************************************************/
/*                             Cell access                                    */
/******************************************************************************/

string& DataTable::operator()(size_t rowIndex, size_t colIndex) noexcept(false)
{
  if (colIndex >= nCol_)
    throw IndexOutOfBoundsException("DataTable::operator(size_t, size_t).", colIndex, 0, nCol_ - 1);
  if (rowIndex >= data_[colIndex].size())
    throw IndexOutOfBoundsException("DataTable::operator(size_t, size_t).", rowIndex, 0, data_[colIndex].size() - 1);
  return data_[colIndex][rowIndex];
}

const string& DataTable::operator()(size_t rowIndex, size_t colIndex) const
noexcept(false)
{
  if (colIndex >= nCol_)
    throw IndexOutOfBoundsException("DataTable::operator(size_t, size_t).", colIndex, 0, nCol_ - 1);
  if (rowIndex >= data_[colIndex].size())
    throw IndexOutOfBoundsException("DataTable::operator(size_t, size_t).", rowIndex, 0, data_[colIndex].size() - 1);
  return data_[colIndex][rowIndex];
}

/******************************************************************************/

string& DataTable::operator()(const string& rowName, const string& colName)
noexcept(false)
{
  if (rowNames_ == NULL)
    throw NoTableRowNamesException("DataTable::operator(const string &, const string &).");
  if (colNames_ == NULL)
    throw NoTableColumnNamesException("DataTable::operator(const string &, const string &).");
  try
  {
    size_t rowIndex = VectorTools::which(*rowNames_, rowName);
    size_t colIndex = VectorTools::which(*colNames_, colName);
    return (*this)(rowIndex, colIndex);
  }
  catch (ElementNotFoundException<string>& ex)
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, const string &).", *ex.getElement());
  }
}

const string& DataTable::operator()(const string& rowName, const string& colName) const
noexcept(false)
{
  if (rowNames_ == NULL)
    throw NoTableRowNamesException("DataTable::operator(const string &, const string &).");
  if (colNames_ == NULL)
    throw NoTableColumnNamesException("DataTable::operator(const string &, const string &).");
  try
  {
    size_t rowIndex = VectorTools::which(*rowNames_, rowName);
    size_t colIndex = VectorTools::which(*colNames_, colName);
    return (*this)(rowIndex, colIndex);
  }
  catch (ElementNotFoundException<string>& ex)
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, const string &).", *ex.getElement());
  }
}

/******************************************************************************/

string& DataTable::operator()(const string& rowName, size_t colIndex)
noexcept(false)
{
  if (rowNames_ == NULL)
    throw NoTableRowNamesException("DataTable::operator(const string &, size_t).");
  if (colIndex >= nCol_)
    throw IndexOutOfBoundsException("DataTable::operator(const string &, size_t).", colIndex, 0, nCol_ - 1);
  try
  {
    size_t rowIndex = VectorTools::which(*rowNames_, rowName);
    return (*this)(rowIndex, colIndex);
  }
  catch (ElementNotFoundException<string>& ex)
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, size_t).", *ex.getElement());
  }
}

const string& DataTable::operator()(const string& rowName, size_t colIndex) const
noexcept(false)
{
  if (rowNames_ == NULL)
    throw NoTableRowNamesException("DataTable::operator(const string &, size_t).");
  if (colIndex >= nCol_)
    throw IndexOutOfBoundsException("DataTable::operator(const string &, size_t).", colIndex, 0, nCol_ - 1);
  try
  {
    size_t rowIndex = VectorTools::which(*rowNames_, rowName);
    return (*this)(rowIndex, colIndex);
  }
  catch (ElementNotFoundException<string>& ex)
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, size_t).", *ex.getElement());
  }
}

/******************************************************************************/

string& DataTable::operator()(size_t rowIndex, const string& colName)
noexcept(false)
{
  if (colNames_ == NULL)
    throw NoTableColumnNamesException("DataTable::operator(size_t, const string &).");
  try
  {
    size_t colIndex = VectorTools::which(*colNames_, colName);
    if (rowIndex >= data_[colIndex].size())
      throw IndexOutOfBoundsException("DataTable::operator(size_t, const string &).", rowIndex, 0, data_[colIndex].size() - 1);
    return (*this)(rowIndex, colIndex);
  }
  catch (ElementNotFoundException<string>& ex)
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, const string &).", *ex.getElement());
  }
}

const string& DataTable::operator()(size_t rowIndex, const string& colName) const
noexcept(false)
{
  if (colNames_ == NULL)
    throw NoTableColumnNamesException("DataTable::operator(size_t, const string &).");
  try
  {
    size_t colIndex = VectorTools::which(*colNames_, colName);
    if (rowIndex >= data_[colIndex].size())
      throw IndexOutOfBoundsException("DataTable::operator(size_t, const string &).", rowIndex, 0, data_[colIndex].size() - 1);
    return (*this)(rowIndex, colIndex);
  }
  catch (ElementNotFoundException<string>& ex)
  {
    throw TableNameNotFoundException("DataTable::operator(const string &, const string &).", *ex.getElement());
  }
}

/******************************************************************************/
/*                             Work with names                                */
/******************************************************************************/

void DataTable::setRowNames(const vector<string>& rowNames)
noexcept(false)
{
  if (!VectorTools::isUnique(rowNames))
  {
    throw DuplicatedTableRowNameException("DataTable::setRowNames(...). Row names must be unique.");
  }
  if (rowNames.size() != nRow_)
    throw DimensionException("DataTable::setRowNames.", rowNames.size(), nRow_);
  else
  {
    if (rowNames_ != NULL)
      delete rowNames_;
    rowNames_ = new vector<string>(rowNames.begin(), rowNames.end());
  }
}

vector<string> DataTable::getRowNames() const noexcept(false)
{
  if (rowNames_ == NULL)
    throw NoTableRowNamesException("DataTable::getRowNames().");
  return *rowNames_;
}

string DataTable::getRowName(size_t index) const noexcept(false)
{
  if (rowNames_ == NULL)
    throw NoTableRowNamesException("DataTable::getRowName(size_t).");
  if (index >= nRow_)
    throw IndexOutOfBoundsException("DataTable::getRowName(size_t).", index, 0, nRow_ - 1);
  return (*rowNames_)[index];
}

/******************************************************************************/

void DataTable::setColumnNames(const vector<string>& colNames)
noexcept(false)
{
  if (!VectorTools::isUnique(colNames))
    throw DuplicatedTableColumnNameException("DataTable::setColumnNames(...). Column names must be unique.");
  if (colNames.size() != nCol_)
    throw DimensionException("DataTable::setColumnNames.", colNames.size(), nCol_);
  else
  {
    if (colNames_ != NULL)
      delete colNames_;
    colNames_ = new vector<string>(colNames.begin(), colNames.end());
  }
}

vector<string> DataTable::getColumnNames() const noexcept(false)
{
  if (colNames_ == NULL)
    throw NoTableColumnNamesException("DataTable::getColumnNames().");
  return *colNames_;
}

string DataTable::getColumnName(size_t index) const noexcept(false)
{
  if (colNames_ == NULL)
    throw NoTableColumnNamesException("DataTable::getColumnName(size_t).");
  if (index >= nCol_)
    throw IndexOutOfBoundsException("DataTable::getColumnName(size_t).", index, 0, nCol_ - 1);
  return (*colNames_)[index];
}

/******************************************************************************/
/*                               Work on columns                              */
/******************************************************************************/

vector<string>& DataTable::getColumn(size_t index)
noexcept(false)
{
  if (index >= nCol_)
    throw IndexOutOfBoundsException("DataTable::getColumn(size_t).", index, 0, nCol_ - 1);
  return data_[index];
}

const vector<string>& DataTable::getColumn(size_t index) const
noexcept(false)
{
  if (index >= nCol_)
    throw IndexOutOfBoundsException("DataTable::getColumn(size_t).", index, 0, nCol_ - 1);
  return data_[index];
}

vector<string>& DataTable::getColumn(const string& colName)
noexcept(false)
{
  if (colNames_ == NULL)
    throw NoTableColumnNamesException("DataTable::getColumn(const string &).");
  try
  {
    size_t colIndex = VectorTools::which(*colNames_, colName);
    return data_[colIndex];
  }
  catch (ElementNotFoundException<string>& ex)
  {
    throw TableColumnNameNotFoundException("DataTable::getColumn(const string &).", colName);
  }
}

const vector<string>& DataTable::getColumn(const string& colName) const
noexcept(false)
{
  if (colNames_ == NULL)
    throw NoTableColumnNamesException("DataTable::getColumn(const string &).");
  try
  {
    size_t colIndex = VectorTools::which(*colNames_, colName);
    return data_[colIndex];
  }
  catch (ElementNotFoundException<string>& ex)
  {
    throw TableColumnNameNotFoundException("DataTable::getColumn(const string &).", colName);
  }
}

bool DataTable::hasColumn(const string& colName) const
{
  if (colNames_ == NULL)
    return false;
  for (size_t i = 0; i < colNames_->size(); i++)
  {
    if ((*colNames_)[i] == colName)
      return true;
  }
  return false;
}

void DataTable::deleteColumn(size_t index)
noexcept(false)
{
  if (index >= nCol_)
    throw IndexOutOfBoundsException("DataTable::deleteColumn(size_t).", index, 0, nCol_ - 1);
  data_.erase(data_.begin() + static_cast<ptrdiff_t>(index));
  if (colNames_)
    colNames_->erase(colNames_->begin() + static_cast<ptrdiff_t>(index));
  nCol_--;
}

void DataTable::deleteColumn(const string& colName)
noexcept(false)
{
  if (!colNames_)
    throw NoTableColumnNamesException("DataTable::deleteColumn(const string &).");
  try
  {
    size_t colIndex = VectorTools::which(*colNames_, colName);
    data_.erase(data_.begin() + static_cast<ptrdiff_t>(colIndex));
    colNames_->erase(colNames_->begin() + static_cast<ptrdiff_t>(colIndex));
    nCol_--;
  }
  catch (ElementNotFoundException<string>& ex)
  {
    throw TableColumnNameNotFoundException("DataTable::deleteColumn(const string &).", colName);
  }
}

void DataTable::addColumn(const vector<string>& newColumn)
noexcept(false)
{
  if (colNames_)
    throw TableColumnNamesException("DataTable::addColumn. Table has column names.");
  if (newColumn.size() != nRow_)
    throw DimensionException("DataTable::addColumn.", newColumn.size(), nRow_);
  data_.push_back(newColumn);
  nCol_++;
}

void DataTable::addColumn(const string& colName, const vector<string>& newColumn)
noexcept(false)
{
  if (!colNames_)
  {
    if (nCol_ == 0)
      colNames_ = new vector<string>();
    else
      throw NoTableColumnNamesException("DataTable::addColumn. Table has column names.");
  }
  if (newColumn.size() != nRow_)
    throw DimensionException("DataTable::addColumn.", newColumn.size(), nRow_);
  if (nCol_ > 0 && find(colNames_->begin(), colNames_->end(), colName) != colNames_->end())
    throw DuplicatedTableColumnNameException("DataTable::addColumn(const string &, const vector<string> &). Column names must be unique.");
  colNames_->push_back(colName);
  data_.push_back(newColumn);
  nCol_++;
}

/******************************************************************************/
/*                               Work on rows                                 */
/******************************************************************************/

vector<string> DataTable::getRow(size_t index) const
noexcept(false)
{
  if (index >= nRow_)
    throw IndexOutOfBoundsException("DataTable::getRow(size_t).", index, 0, nRow_ - 1);
  vector<string> row;
  for (size_t i = 0; i < nCol_; i++)
  {
    row.push_back(data_[i][index]);
  }
  return row;
}

vector<string> DataTable::getRow(const string& rowName) const
noexcept(false)
{
  if (!rowNames_)
    throw NoTableRowNamesException("DataTable::getRow(const string &).");
  try
  {
    size_t rowIndex = VectorTools::which(*rowNames_, rowName);
    vector<string> row;
    for (size_t i = 0; i < nCol_; i++)
    {
      row.push_back(data_[i][rowIndex]);
    }
    return row;
  }
  catch (ElementNotFoundException<string>& ex)
  {
    throw TableRowNameNotFoundException("DataTable::getRow(const string &).", rowName);
  }
}

bool DataTable::hasRow(const string& rowName) const
{
  if (rowNames_ == NULL)
    return false;
  for (size_t i = 0; i < rowNames_->size(); i++)
  {
    if ((*rowNames_)[i] == rowName)
      return true;
  }
  return false;
}

void DataTable::deleteRow(size_t index)
noexcept(false)
{
  for (size_t j = 0; j < nCol_; j++)
  {
    vector<string>* column = &data_[j];
    if (index >= column->size())
      throw IndexOutOfBoundsException("DataTable::deleteRow(size_t).", index, 0, column->size() - 1);
    column->erase(column->begin() + static_cast<ptrdiff_t>(index));
  }
  if (rowNames_)
    rowNames_->erase(rowNames_->begin() + static_cast<ptrdiff_t>(index));
  nRow_--;
}

void DataTable::deleteRow(const string& rowName)
noexcept(false)
{
  if (!rowNames_)
    throw NoTableRowNamesException("DataTable::deleteRow(const string &).");
  try
  {
    size_t rowIndex = VectorTools::which(*rowNames_, rowName);
    for (size_t j = 0; j < nCol_; j++)
    {
      vector<string>* column = &data_[j];
      column->erase(column->begin() + static_cast<ptrdiff_t>(rowIndex));
    }
    rowNames_->erase(rowNames_->begin() + static_cast<ptrdiff_t>(rowIndex));
    nRow_--;
  }
  catch (ElementNotFoundException<string>& ex)
  {
    throw TableRowNameNotFoundException("DataTable::deleteRow(const string &).", rowName);
  }
}

void DataTable::addRow(const vector<string>& newRow)
noexcept(false)
{
  if (rowNames_)
    throw TableRowNamesException("DataTable::addRow. Table has row names.");
  if (newRow.size() != nCol_)
    throw DimensionException("DataTable::addRow.", newRow.size(), nCol_);
  for (size_t j = 0; j < nCol_; j++)
  {
    data_[j].push_back(newRow[j]);
  }
  nRow_++;
}

void DataTable::addRow(const string& rowName, const vector<string>& newRow)
noexcept(false)
{
  if (!rowNames_)
  {
    if (nRow_ == 0)
      rowNames_ = new vector<string>();
    else
      throw NoTableRowNamesException("DataTable::addRow. Table has row names.");
  }
  if (newRow.size() != nCol_)
    throw DimensionException("DataTable::addRow.", newRow.size(), nCol_);
  if (nRow_ > 0 && find(rowNames_->begin(), rowNames_->end(), rowName) != rowNames_->end())
    throw DuplicatedTableRowNameException("DataTable::addRow(const string &, const vector<string> &). Row names must be unique.");
  rowNames_->push_back(rowName);
  for (size_t j = 0; j < nCol_; j++)
  {
    data_[j].push_back(newRow[j]);
  }
  nRow_++;
}

/******************************************************************************/
/*                               Read from a CSV file                         */
/******************************************************************************/

DataTable* DataTable::read(istream& in, const string& sep, bool header, int rowNames)
noexcept(false)
{
  string firstLine  = FileTools::getNextLine(in);
  StringTokenizer st1(firstLine, sep, false, true);
  vector<string> row1(st1.getTokens().begin(), st1.getTokens().end());
  string secondLine = FileTools::getNextLine(in);
  StringTokenizer st2(secondLine, sep, false, true);
  vector<string> row2(st2.getTokens().begin(), st2.getTokens().end());
  size_t nCol = row1.size();
  bool hasRowNames;
  DataTable* dt;
  if (row1.size() == row2.size())
  {
    dt = new DataTable(nCol);
    if (header)
    { // Use first line as header.
      dt->setColumnNames(row1);
    }
    else
    {
      dt->addRow(row1);
    }
    dt->addRow(row2);
    hasRowNames = false;
  }
  else if (row1.size() == row2.size() - 1)
  {
    dt = new DataTable(nCol);
    dt->setColumnNames(row1);
    string rowName = *row2.begin();
    dt->addRow(rowName, vector<string>(row2.begin() + 1, row2.end()));
    hasRowNames = true;
  }
  else
    throw DimensionException("DataTable::read(...). Row 2 has not the correct number of columns.", row2.size(), nCol);

  // Now read each line:
  string line = FileTools::getNextLine(in);
  while (!TextTools::isEmpty(line))
  {
    StringTokenizer st(line, sep, false, true);
    if (hasRowNames)
    {
      string rowName = *st.getTokens().begin();
      vector<string> row(st.getTokens().begin() + 1, st.getTokens().end());
      dt->addRow(rowName, row);
    }
    else
    {
      vector<string> row(st.getTokens().begin(), st.getTokens().end());
      dt->addRow(row);
    }
    line = FileTools::getNextLine(in);
  }

  // Row names:
  if (rowNames > -1)
  {
    if (static_cast<size_t>(rowNames) >= nCol)
      throw IndexOutOfBoundsException("DataTable::read(...). Invalid column specified for row names.", static_cast<size_t>(rowNames), 0, nCol - 1);
    vector<string> col = dt->getColumn(static_cast<size_t>(rowNames));
    dt->setRowNames(col);
    dt->deleteColumn(static_cast<size_t>(rowNames));
  }

  return dt;
}

/******************************************************************************/

void DataTable::write(const DataTable& data, ostream& out, const string& sep, bool alignHeaders)
{
  size_t n = data.getNumberOfColumns();
  if (n == 0)
    return;
  if (data.hasColumnNames())
  { // Write header
    vector<string> colNames = data.getColumnNames();
    if (alignHeaders && data.hasRowNames())
      out << sep;
    out << colNames[0];
    for (size_t i = 1; i < n; i++)
    {
      out << sep << colNames[i];
    }
    out << endl;
  }
  // Now write each row:
  for (size_t i = 0; i < data.getNumberOfRows(); i++)
  {
    if (data.hasRowNames())
      out << data.getRowName(i) << sep;
    out << data(i, 0);
    for (size_t j = 1; j < n; j++)
    {
      out << sep << data(i, j);
    }
    out << endl;
  }
}

void DataTable::write(const DataTable& data, bpp::OutputStream& out, const string& sep, bool alignHeaders)
{
  size_t n = data.getNumberOfColumns();
  if (n == 0)
    return;
  if (data.hasColumnNames())
  { // Write header
    vector<string> colNames = data.getColumnNames();
    if (alignHeaders && data.hasRowNames())
      out << sep;
    out << colNames[0];
    for (size_t i = 1; i < n; i++)
    {
      out << sep << colNames[i];
    }
    out.endLine();
  }
  // Now write each row:
  for (size_t i = 0; i < data.getNumberOfRows(); i++)
  {
    if (data.hasRowNames())
      out << data.getRowName(i) << sep;
    out << data(i, 0);
    for (size_t j = 1; j < n; j++)
    {
      out << sep << data(i, j);
    }
    out.endLine();
  }
}

/******************************************************************************/

