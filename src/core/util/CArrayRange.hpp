#pragma once


/**
 * A wrapper arround C-style array to allow range iteration:
 *  int *array = new[size];
 *  ...
 *  for (int elem: CArrayRange<int>(array, size)) {
 *    // do something
 *  }
 */
template <typename DataType>
class CArrayRange {
public:
 class iterator {
 public:
   iterator(DataType * ptr): ptr(ptr){}
   iterator operator++() { ++ptr; return *this; }
   bool operator!=(const iterator & other) const { return ptr != other.ptr; }
   const DataType& operator*() const { return *ptr; }
 private:
   DataType* ptr;
 };
private:
 DataType *val;
 unsigned len;
public:
 CArrayRange(DataType *ptr, size_t len): val(ptr), len(len) {}
 iterator begin() const { return iterator(val); }
 iterator end() const { return iterator(val + len); }
 unsigned int size() const { return len; }
};

