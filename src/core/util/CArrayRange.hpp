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
   iterator(DataType * ptr): _ptr(ptr){}
   iterator operator++() { ++_ptr; return *this; }
   bool operator!=(const iterator & other) const { return _ptr != other._ptr; }
   const DataType& operator*() const { return *_ptr; }
 private:
   DataType* _ptr;
 };
private:
 DataType *_val;
 size_t _len;
public:
 CArrayRange(DataType *ptr, size_t len): _val(ptr), _len(len) {}
 iterator begin() const { return iterator(_val); }
 iterator end() const { return iterator(_val + _len); }
 size_t size() const { return _len; }
};

