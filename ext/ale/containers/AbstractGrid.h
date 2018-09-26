// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 19/10/16.
//

#ifndef GRID_H
#define GRID_H

//Includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <cassert>
#include <fstream>
#include <sstream>


template<class T>
std::size_t getStreamElementSize(const T& o, const std::ostream& os) {
  /// Returns length in stream of a given element.
  std::ostringstream temp;
  //temp.precision(2);
  auto old_precision = os.precision();
  //os.precision(old_precision);
  temp << std::setprecision(old_precision) << o;
  return temp.str().size();
}


template <typename T, typename rowKey, typename colKey>
class AbstractGrid {
  /*!
   * \class AbstractGrid
   * \brief Abstract grid, a 2D vector (nrow_, ncol_) of T objects.
   * \details AbstractGrid is a template class whith functions to define. These functions gives the behaviour of the 2D
   * grid as the kinf of indexes to use, how to find a value according to the index etc.
   * See Grid for a simple definition of a Matrix using AbstractGrid.
   *
   */
private:
/****************
 * Private getters
 */
  inline T& _get(const std::size_t x, const std::size_t y);

/****************
 * Private setters
 */
  inline const T& _get(const std::size_t x, const std::size_t y) const;

  // Special getters to vectors
  inline std::vector<T> _getRow(const std::size_t x) const;
  inline std::vector<T> _getCol(const std::size_t y) const;

  // Setters
  void _set(const std::size_t x, const std::size_t y, const T& value);
  void _removeRow(const std::size_t i);
  void _removeCol(const std::size_t i);
  void _addRows_(const AbstractGrid<T, rowKey, colKey> &x, const std::vector<rowKey> &keys);
  void _addCols_(const AbstractGrid<T, rowKey, colKey> &y, const std::vector<colKey> &keys);
  void _addRow(const std::vector<T> &x);
  void _addCol(const std::vector<T> &y);

protected:
  /// Number of rows.
  std::size_t nrow_;
  /// Number of columns.
  std::size_t ncol_;
  /// 2D vector which contains data.
  std::vector<std::vector<T> > data_;

  /// Returns the unsigned integer associated to the key.
  virtual std::size_t _rowKeyToIndex(const rowKey& key) const = 0;
  /// Returns the unsigned integer associated to the key.
  virtual std::size_t _colKeyToIndex(const colKey& key) const = 0;
  /// Returns true if the row index exists, false if not.
  virtual bool _isRowKey(const rowKey& key) const = 0;
  /// Returns true if the col index exists, false if not.
  virtual bool _isColKey(const colKey& key) const = 0;
  /// Get rowKeys
  virtual std::vector<rowKey> _getRowKeys() const = 0;
  /// Get colKeys
  virtual std::vector<colKey> _getColKeys() const = 0;
  /// Add a row index.
  virtual void _addRowKey(const rowKey &key, const std::size_t &index) = 0;
  /// Add a col index.
  virtual void _addColKey(const colKey &key, const std::size_t &index) = 0;
  /// Erase a row key.
  virtual void _removeRowKey(const rowKey& key) = 0;
  /// Erase a col key.
  virtual void _removeColKey(const colKey& key) = 0;
  /// Erase all row keys.
  virtual void _eraseAllRowKeys() = 0;
  /// Erase all col keys.
  virtual void _eraseAllColKeys() = 0;

public:
/****************
 * Constructors
 */
  AbstractGrid();

  AbstractGrid(const std::size_t nrow, const std::size_t ncol);

  AbstractGrid(const std::size_t nrow, const std::size_t ncol, const T& model);

  explicit AbstractGrid(const std::vector<std::vector<T> >& grid);

  AbstractGrid(const AbstractGrid<T, rowKey, colKey>& grid);

  AbstractGrid(AbstractGrid<T, rowKey, colKey>&& grid) noexcept;

/****************
 * Destructor
 */
  virtual ~AbstractGrid() { this->clear(); }

/****************
 * Getters
 */
  virtual T& get(const rowKey& i, const colKey& j) { return _get(_rowKeyToIndex(i), _colKeyToIndex(j)); }

  T& operator()(const rowKey& i, const colKey& j) { return get(i, j); }

  std::vector<T>& operator[](const rowKey& i) { return data_[_rowKeyToIndex(i)]; }

  virtual const T& get(const rowKey& i, const colKey& j) const { return _get(_rowKeyToIndex(i), _colKeyToIndex(j)); }

  const T& operator()(const rowKey& i, const colKey& j) const { return get(i, j); }

  const std::vector<T>& operator[](const rowKey& i) const { return data_.at(_rowKeyToIndex(i)); }

  std::vector<T> getRow(const rowKey& i) const { return _getRow(_rowKeyToIndex(i)); }

  std::vector<T> getCol(const colKey& j) const { return _getCol(_colKeyToIndex(j)); }

  std::size_t ncol() const { return ncol_; }

  std::size_t nrow() const { return nrow_; }

  std::vector<std::vector<T>> to2DVector() const { return data_; }

  std::vector<std::vector<T>> toStl() const { return data_; }

  std::vector<rowKey> getRowIndexes() const { return _getRowKeys(); }

  std::vector<colKey> getColIndexes() const { return _getColKeys(); }

  template<typename Comp>
  bool equalTo(const AbstractGrid<T, rowKey, colKey>& other, const Comp& comp) const {
    /// Compare the current AbstractGrid to an other.
    if(!(this->nrow_ == other.nrow_ and this->ncol_ == other.ncol_)){
      return false;
    }

    auto this_rowKeys = getRowIndexes();
    auto this_colKeys = getColIndexes();
    for(auto& this_rowKey: this_rowKeys){
      if(not other.hasRowIndex(this_rowKey))
        return false;
      for(auto& this_colKey: this_colKeys){
        if(not other.hasColIndex(this_colKey))
          return false;
        const auto& left_elem = this->operator()(this_rowKey, this_colKey);
        const auto& right_elem = other(this_rowKey, this_colKey);
        if(not comp(left_elem, right_elem))
          return false;
      }
    }

    return true;
  }

  void print() const;

  std::pair<rowKey, colKey> min() const;

  std::pair<rowKey, colKey> max() const;

  bool hasRowIndex(const rowKey& key) const { return _isRowKey(key); }

  bool hasColIndex(const colKey& key) const { return _isColKey(key); }

  void csv(const std::string& filename = "savefile.csv") const;

  template<class CLASS>
  void write(std::ostream& os, const CLASS& e, const std::size_t max_length = 5) const {
    auto size = getStreamElementSize(e, os);
    os << e;
    while ((size++) < max_length) os << ' ';
  }

  std::size_t maxStreamSize(const std::ostream& os) const;

  /// Print a Table into an output stream.
  void write(std::ostream& os
             , const std::string& sep = "\t"
             , const bool rownames = true
             , const bool colnames = true
             , const bool align = true) const;

/****************
 * Setters
 */
  AbstractGrid<T, rowKey, colKey>& operator=(const AbstractGrid<T, rowKey, colKey>& right);

  AbstractGrid<T, rowKey, colKey>& operator=(AbstractGrid<T, rowKey, colKey>&& right) noexcept;

  AbstractGrid<T, rowKey, colKey>& operator=(const std::vector<std::vector<T>>& right);

  AbstractGrid<T, rowKey, colKey>& operator=(const std::vector<T>& right);

  void set(const rowKey& i, const colKey& j, const T& value) { _set(_rowKeyToIndex(i), _colKeyToIndex(j), value); }

  void addRow(const std::vector<T> &x, const rowKey& key);

  void addCol(const std::vector<T> &y, const colKey& key);

  void addRows(const AbstractGrid<T, rowKey, colKey> &x, const std::vector<rowKey> &keys);

  void addCols(const AbstractGrid<T, rowKey, colKey> &y, const std::vector<colKey> &keys);

  void removeRow(const rowKey& row) { _removeRow(_rowKeyToIndex(row)); _removeRowKey(row); }

  void removeCol(const colKey& col) { _removeCol(_colKeyToIndex(col)); _removeColKey(col); }

  void clear();

  friend bool operator!=(const AbstractGrid<T, rowKey, colKey> &lhs
                         , const AbstractGrid<T, rowKey, colKey> &lhr
  ){
    return !(lhs == lhr);
  }

  friend bool operator==(const AbstractGrid<T, rowKey, colKey> &lhs
                         , const AbstractGrid<T, rowKey, colKey> &lhr
  ){
    return lhs.equalTo(lhr, [](const T& a, const T& b){ return a == b; });
  }

  friend std::ostream& operator<<(std::ostream& os
                                  , const AbstractGrid<T, rowKey, colKey>& grid
  ) {
    grid.write(os, "\t", true, true); return os;
  }
};

/****************
 * Implementations
 */
template <typename T, typename rowKey, typename colKey>
AbstractGrid<T, rowKey, colKey>::AbstractGrid() : nrow_(0), ncol_(0), data_(std::vector<std::vector<T>>(0, std::vector<T>(0))) {
}

template <typename T, typename rowKey, typename colKey>
AbstractGrid<T, rowKey, colKey>::AbstractGrid(const std::size_t nrow, const std::size_t ncol) : nrow_(nrow), ncol_(ncol) {
  data_ = std::vector<std::vector<T> >(nrow, std::vector<T>(ncol));
}

template <typename T, typename rowKey, typename colKey>
AbstractGrid<T, rowKey, colKey>::AbstractGrid(const std::size_t nrow, const std::size_t ncol, const T& model) : nrow_(nrow), ncol_(ncol) {
  data_ = std::vector<std::vector<T> >(nrow, std::vector<T>(ncol, model));
}

template <typename T, typename rowKey, typename colKey>
AbstractGrid<T, rowKey, colKey>::AbstractGrid(const std::vector<std::vector<T> >& grid) {
  nrow_ = grid.size();
  ncol_ = grid[0].size();
  data_ = grid;
}

template <typename T, typename rowKey, typename colKey>
AbstractGrid<T, rowKey, colKey>::AbstractGrid(const AbstractGrid<T, rowKey, colKey>& grid) {
  nrow_ = grid.nrow_;
  ncol_ = grid.ncol_;
  data_ = grid.data_;
}

template <typename T, typename rowKey, typename colKey>
AbstractGrid<T, rowKey, colKey>::AbstractGrid(AbstractGrid<T, rowKey, colKey>&& grid) noexcept {
  nrow_ = std::exchange(grid.nrow_, 0);
  ncol_ = std::exchange(grid.ncol_, 0);
  data_ = std::move(grid.data_);
}

template <typename T, typename rowKey, typename colKey>
void AbstractGrid<T, rowKey, colKey>::_set(const std::size_t x, const std::size_t y, const T& value){
  data_[x][y] = value;
}

template <typename T, typename rowKey, typename colKey>
void AbstractGrid<T, rowKey, colKey>::print() const {
  /// Obsolete method to print the 2D grid in terminal. See write(...) method.
  if(nrow_ == 0)
    std::cout<<"The grid is empty"<<std::endl;
  else {
    for (std::size_t i = 0 ; i < nrow_; ++i){
      for (std::size_t j = 0 ; j < ncol_; ++j) {
        std::cout << data_[i][j] << " ";
      }
      std::cout << std::endl;
    }
  }
}

//non-const getters
template <typename T, typename rowKey, typename colKey>
T& AbstractGrid<T, rowKey, colKey>::_get(const std::size_t x, const std::size_t y){
  if(x >= nrow_ or y >= ncol_) {
    std::cerr << "Error: bad index:" << std::endl;
    if(x >= nrow_) std::cerr << "\tLine index = " << x << " is not valid." << std::endl;
    else std::cerr << "\tColumn index = " << y << " is not valid." << std::endl;
    exit(EXIT_FAILURE);
  }
  return data_[x][y];
}

//const getters
template <typename T, typename rowKey, typename colKey>
const T& AbstractGrid<T, rowKey, colKey>::_get(const std::size_t x, const std::size_t y) const{
  if(x >= nrow_ or y >= ncol_) {
    std::cerr << "Error: bad index:" << std::endl;
    if(x >= nrow_) std::cerr << "\tLine index = " << x << " is not valid." << std::endl;
    else std::cerr << "\tColumn index = " << y << " is not valid." << std::endl;
    exit(EXIT_FAILURE);
  }
  return data_[x][y];
}

template<typename T, typename rowKey, typename colKey>
std::vector<T> AbstractGrid<T, rowKey, colKey>::_getRow(const std::size_t x) const {
  if(x >= nrow_) {
    std::cerr << "Error: bad index:" << std::endl;
    std::cerr << "\tLine index = " << x << " is not valid." << std::endl;
    exit(EXIT_FAILURE);
  }
  return data_[x];
}

template<typename T, typename rowKey, typename colKey>
std::vector<T> AbstractGrid<T, rowKey, colKey>::_getCol(const std::size_t y) const {
  if(y >= ncol_) {
    std::cerr << "Error: bad index." << std::endl;
    std::cerr << "\tColumn index = " << y << " is not valid." << std::endl;
    exit(EXIT_FAILURE);
  }
  std::vector<T> ret;
  ret.reserve(nrow_);
  for(std::size_t i = 0 ; i < nrow_ ; ++i){
    ret.push_back(data_[i][y]);
  }
  return ret;
}

template<typename T, typename rowKey, typename colKey>
void AbstractGrid<T, rowKey, colKey>::_addRow(const std::vector<T> &x) {
  if((x.size() == this->ncol_) or (this->nrow_ == 0 and this->ncol_ == 0)) {
    if(this->ncol_ == 0)
      this->ncol_ = x.size();
    this->data_.push_back(x);
    this->nrow_+=1;
  }
  else {
    std::cerr << "AbstractGrid<T>::_addRow( ... ) : bad vector size with " << x.size() << " elements, " << ncol_
              << " expected." << std::endl;
    exit(EXIT_FAILURE);
  }
}

template<typename T, typename rowKey, typename colKey>
void AbstractGrid<T, rowKey, colKey>::_addCol(const std::vector<T> &y) {
  if(y.size() == this->nrow_ or (this->nrow_ == 0 and this->ncol_ == 0)){
    if(this->nrow_ == 0) {
      this->nrow_ = y.size();
      //data_.reserve(this->nrow_);
      data_ = std::vector<std::vector<T>>(y.size());//, std::vector<T>());
    }
    for(std::size_t i = 0; i < this->nrow_ ; ++i){
      this->data_[i].emplace_back(y.at(i));
    }
    this->ncol_+=1;
  }
  else {
    std::cerr << "AbstractGrid<T>::_addCol( ... ) : bad vector size with " << y.size() << " elements, " << nrow_ << " expected." << std::endl;
    exit(EXIT_FAILURE);
  }
}

template<typename T, typename rowKey, typename colKey>
inline void
AbstractGrid<T, rowKey, colKey>::_addRows_(const AbstractGrid<T, rowKey, colKey> &x, const std::vector<rowKey> &keys) {
  for(std::size_t i = 0; i < x.nrow_; ++i) {
    this->_addRow(x[i]);
    this->_addRowKey(keys.at(i));
  }
}

template<typename T, typename rowKey, typename colKey>
inline void
AbstractGrid<T, rowKey, colKey>::_addCols_(const AbstractGrid<T, rowKey, colKey> &y, const std::vector<colKey> &keys) {
  if(y.nrow_ == nrow_ or (nrow_ == 0 and ncol_ == 0)) {
    for (std::size_t i = 0; i < y.nrow_; ++i) {
      auto &grid_line = data_[i];
      auto &inser_line = y[i];
      grid_line.insert(grid_line.end(),
                       inser_line.begin(), inser_line.end());
      this->_addColKey(keys.at(i));
    }
    ncol_+=y.ncol_;
  }
  else {
    std::cerr << "AbstractGrid<T, rowKey, colKey>::_addCols_( ... ) : bad vector size with " << y.nrow_ << " elements, " << nrow_ << " expected." << std::endl;
    exit(EXIT_FAILURE);
  }
}

template<typename T, typename rowKey, typename colKey>
void
AbstractGrid<T, rowKey, colKey>::addRows(const AbstractGrid<T, rowKey, colKey> &x, const std::vector<rowKey> &keys) {
  /// Add rows of an AbstractGrid.
  if(&x != this){// if the AbstractGrid is itself
    _addRows_(x, keys);
  } else {
    AbstractGrid<T, rowKey, colKey> copy(x);
    this->_addRows_(copy, keys);
  }
}

template<typename T, typename rowKey, typename colKey>
void
AbstractGrid<T, rowKey, colKey>::addCols(const AbstractGrid<T, rowKey, colKey> &y, const std::vector<colKey> &keys) {
  /// Add cols of an AbstractGrid.
  if(&y != this){// if the AbstractGrid is itself
    _addRows_(y, keys);
  } else {
    AbstractGrid<T, rowKey, colKey> copy(y);
    _addCols_(copy, keys);
  }
}

template <typename T, typename rowKey, typename colKey>
AbstractGrid<T, rowKey, colKey>& AbstractGrid<T, rowKey, colKey>::operator=(const AbstractGrid<T, rowKey, colKey>& right){
  if(&right != this) {
    this->nrow_ = right.nrow_;
    this->ncol_ = right.ncol_;
    this->data_ = right.data_;
  }
  return (*this);
}

template <typename T, typename rowKey, typename colKey>
AbstractGrid<T, rowKey, colKey>& AbstractGrid<T, rowKey, colKey>::operator=(AbstractGrid<T, rowKey, colKey>&& right) noexcept {
  if(&right != this) {
    this->nrow_ = std::exchange(right.nrow_, 0);
    this->ncol_ = std::exchange(right.ncol_, 0);
    this->data_ = std::move(right.data_);
  }
  return (*this);
}

template <typename T, typename rowKey, typename colKey>
AbstractGrid<T, rowKey, colKey>& AbstractGrid<T, rowKey, colKey>::operator=(const std::vector<std::vector<T>>& right){
  /// Assign a 2d std::vector to build an AbstractGrid. Use with precautions.
  for(std::size_t i = 0; i < right.size(); ++i){
    for(std::size_t j = 0; j < right.at(i).size(); ++j){
      _addColKey(colKey(i), i);
    }
    _addRowKey(rowKey(i), i);
  }
  this->nrow_ = right.size();
  if(right.size() != 0)
    this->ncol_ = right[0].size();
  else
    this->ncol_ = 0;
  this->data_ = right;
  return (*this);
}

template <typename T, typename rowKey, typename colKey>
AbstractGrid<T, rowKey, colKey>& AbstractGrid<T, rowKey, colKey>::operator=(const std::vector<T>& right){
  /// Assign a std::vector as row to build an AbstractGrid(1, vector.size()). Use with precautions.
  this->clear();
  _addRowKey(rowKey(0), 0);
  for(std::size_t i = 0; i < right.size(); ++i){
    _addColKey(colKey(i), i);
  }
  this->_addRow(right);
  return (*this);
}

template <typename T, typename rowKey, typename colKey>
void AbstractGrid<T, rowKey, colKey>::_removeRow(const std::size_t i) {
  /// Remove a row but not its key in AbstractGrid.
  if(this->nrow_ > 0) {
    this->data_.erase(data_.begin() + i);
    this->nrow_ -= 1;
  }
  if(this->nrow_ == 0) this->ncol_ = 0;
  assert(nrow_ >= 0);
}

template <typename T, typename rowKey, typename colKey>
void AbstractGrid<T, rowKey, colKey>::_removeCol(const std::size_t i) {
  /// Remove a column but not its key in AbstractGrid.
  if(this->ncol_ > 1) {
    for(auto itline = this->data_.begin(); itline != this->data_.end(); itline++){
      itline->erase(itline->begin() + i);
      if(itline->empty()){
        this->_removeRow(itline - data_.begin());
      }
    }
    this->ncol_ -= 1;
  }
  else if(this->ncol_ == 1){
    this->clear();
  }
  assert(ncol_ >= 0);
}

template<typename T, typename rowKey, typename colKey>
void AbstractGrid<T, rowKey, colKey>::clear() {
  /// Delete all elements of the container.
  for(auto& e: data_)
    e.clear();

  data_.clear();
  this->nrow_ = 0;
  this->ncol_ = 0;
}

template <typename T, typename rowKey, typename colKey>
void AbstractGrid<T, rowKey, colKey>::addRow(const std::vector<T> &x, const rowKey &key) {
  /// Add a row.
  if(_isRowKey(key))
    removeRow(key);
  _addRowKey(key, nrow_);
  _addRow(x);
}

template <typename T, typename rowKey, typename colKey>
void AbstractGrid<T, rowKey, colKey>::addCol(const std::vector<T> &x, const colKey &key) {
  /// Add column.
  if(_isColKey(key))
    removeCol(key);
  _addColKey(key, ncol_);
  _addCol(x);
}

template <typename T, typename rowKey, typename colKey>
std::pair<rowKey, colKey> AbstractGrid<T, rowKey, colKey>::min() const {
  /// Returns the position of the minimal value.
  auto rowKeys = getRowIndexes();
  auto colKeys = getColIndexes();
  rowKey x_min = rowKeys.front();
  colKey y_min = colKeys.front();
  T min = this->get(x_min, y_min);
  for(auto& x: rowKeys){
    for(auto& y: colKeys){
      if(this->get(x, y) < min){
        x_min = x;
        y_min = y;
        min = this->_get(x_min, y_min);
      }
    }
  }
  return std::pair<rowKey, colKey>(x_min, y_min);
}

template <typename T, typename rowKey, typename colKey>
std::pair<rowKey, colKey> AbstractGrid<T, rowKey, colKey>::max() const {
  /// Returns the position of the maximal value.
  auto rowKeys = getRowIndexes();
  auto colKeys = getColIndexes();
  rowKey& x_max = rowKeys.front();
  colKey& y_max = colKeys.front();
  T max = this->get(x_max, y_max);
  for(auto& x: rowKeys){
    for(auto& y: colKeys){
      if(this->get(x, y) > max){
        x_max = x;
        y_max = y;
        max = this->_get(x_max, y_max);
      }
    }
  }
  return std::pair<rowKey, colKey>(x_max, y_max);
}

template <typename T, typename rowKey, typename colKey>
void AbstractGrid<T, rowKey, colKey>::csv(const std::string &filename) const {
  /// Save the grid as a CSV file.
  std::ofstream file(filename);
  this->write(file, ",", true, true);
  file.close();
}

template<typename T, typename rowKey, typename colKey>
std::size_t AbstractGrid<T, rowKey, colKey>::maxStreamSize(const std::ostream& os) const {
  /// Returns the maximum length of an element in stream.
  std::size_t max_length = 0;
  std::size_t length;
  auto rownames = this->getRowIndexes();
  auto colnames = this->getColIndexes();
  for(std::size_t i = 0; i < this->nrow_; i++){
    length = getStreamElementSize(rownames.at(i), os);
    if(max_length < length)
      max_length = length;

    for(std::size_t j = 0; j < this->ncol_; j++){
      if(i == 0){
        length = getStreamElementSize(colnames.at(i), os);
        if(max_length < length)
          max_length = length;
      }
      auto e = this->get(rownames.at(i),
                         colnames.at(j));
      length = getStreamElementSize(e, os);
      if(max_length < length)
        max_length = length;
    }
  }

  return max_length;
}

template<typename T, typename rowKey, typename colKey>
void AbstractGrid<T, rowKey, colKey>::write(std::ostream &os, const std::string &sep, const bool rownames, const bool colnames,
                         const bool align) const {
  /// Print content and indexes (default) in a given output stream.
  std::size_t max_length = (align) ? maxStreamSize(os) : 0;

  auto rowIndexes = this->getRowIndexes();

  auto colIndexes = this->getColIndexes();

  if(colnames) {
    this->write(os, "", max_length);
    os << sep;
    for (auto colnames_it = colIndexes.begin(); colnames_it != colIndexes.end(); colnames_it++) {
      this->write(os, *colnames_it, max_length);
      if(colnames_it != (colIndexes.end()-1)) os << sep;
    }
    os << std::endl;
  }

  for(auto rowindexes_it = rowIndexes.begin(); rowindexes_it!= rowIndexes.end(); rowindexes_it++){
    if(rownames) {
      this->write(os, *rowindexes_it, max_length);
      os << sep;
    }
    for(auto colindexes_it = colIndexes.begin(); colindexes_it!= colIndexes.end(); colindexes_it++) {
      auto& e = (*this)(*rowindexes_it, *colindexes_it);
      this->write(os, e, max_length);
      if(colindexes_it != (colIndexes.end()-1)) os << sep;
    }
    os << std::endl;
  }
}

#endif //GRID_H

