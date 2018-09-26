//
// File: Range.h
// Created by: Julien Dutheil
// Created on: Mon Nov 21 15:52 2011
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

#ifndef _RANGE_H_
#define _RANGE_H_

#include "../Text/TextTools.h"
#include "../Clonable.h"

//From the STL:
#include <string>
#include <set>
#include <algorithm>
#include <iostream>
#include <cstddef>

namespace bpp {

/**
 * @brief The Range class, defining an interval.
 *
 * Methods are provided for extending the range, get union and intersection.
 */
template<class T> class Range:
  public virtual Clonable
{
  private:
    T begin_;
    T end_;

  public:
    /**
     * @brief Creates a new interval.
     *
     * If a > b, then the positions are swapped.
     * If a == b, the interval is considered empty.
     * Coordinates are 0-based and of type [a, b[,
     * so that the length of the interval is computed as
     * b - a.
     *
     * @param a First position
     * @param b Second position
     */
    Range(const T& a = 0, const T& b = 0):
      begin_(std::min(a, b)),
      end_(std::max(a, b))
    {}

    Range(const Range<T>& range): begin_(range.begin_), end_(range.end_) {}
    
    Range<T>& operator=(const Range<T>& range) {
      begin_ = range.begin_;
      end_   = range.end_;
      return *this;
    }

    Range<T>* clone() const { return new Range<T>(*this); }

    virtual ~Range() {}

  public:
    bool operator==(const Range<T>& r) const {
      return begin_ == r.begin_ && end_ == r.end_;
    }
    bool operator!=(const Range<T>& r) const {
      return begin_ != r.begin_ || end_ != r.end_;
    }
    bool operator<(const Range<T>& r) const {
      return begin_ < r.begin_ || end_ < r.end_;
    }
    virtual Range& operator+=(const T& val) {
      begin_ += val;
      end_ += val;
      return *this;
    }
    virtual Range operator+(const T& val) {
      return Range<T>(*this) += val;
    }
    virtual Range& operator-=(const T& val) {
      begin_ -= val;
      end_ -= val;
      return *this;
    }
    virtual Range operator-(const T& val) {
      return Range<T>(*this) -= val;
    }

    T begin() const { return begin_; }
    
    T end() const { return end_; }
    
    T length() const { return end_ - begin_; }

    /**
     * @param r Range to compare with.
     * @return True if the two intervals overlap.
     */
    bool overlap(const Range& r) const
    {
      return (r.begin_ < end_ && r.end_ > begin_);
    }

    /**
     * @param r Range to compare with.
     * @return True if the two intervals are contiguous (i.e. the two intervals
     * are adjacent and share one bound).
     */
    bool isContiguous(const Range& r) const
    {
      return (r.begin_ == end_ || r.end_ == begin_);
    }

    /**
     * @param r Range to compare with.
     * @return True if the given interval is included in the instanciated one.
     */
    bool contains(const Range& r) const
    {
      return (r.begin_ >= begin_ && r.end_ <= end_);
    }

    /**
     * @brief Expand the current interval with the given one.
     *
     * If the two intervals do not overlap, then the interval is not modified.
     * @param r input interval.
     */
    void expandWith(const Range& r)
    {
      if (r.begin_ < begin_ && r.end_ >= begin_)
        begin_ = r.begin_;
      if (r.end_ > end_ && r.begin_ <= end_)
        end_ = r.end_;
    }

    /**
     * @brief Restrict the current interval to the intersection with the given one.
     *
     * If the two intervals do not overlap, then the interval is set to empty.
     * @param r input interval.
     */
    void sliceWith(const Range& r)
    {
      if (!overlap(r)) {
        begin_ = 0;
        end_   = 0;
      } else {
        if (r.begin_ > begin_ && r.begin_ <= end_)
          begin_ = r.begin_;
        if (r.end_ < end_ && r.end_ >= begin_)
          end_ = r.end_;
      }
    }

    /**
     * @return True if then begining position equals the ending one.
     */
    bool isEmpty() const { return begin_ == end_; }

    /**
     * @return A string describing the range.
     */
    std::string toString() const {
      return ("[" + TextTools::toString(begin_) + "," + TextTools::toString(end_) + "[");
    }

};

/**
 * @brief Interface discribing a collection of Range objects.
 */
template<class T> class RangeCollection {
  public:
    virtual ~RangeCollection() {}
    /**
     * @brief Add a new range to the collection.
     *
     * @param r The range to add to the collection.
     */
    virtual void addRange(const Range<T>& r) = 0;

    /**
     * @brief Get the intersection with a given range.
     *
     * The new multirange is the union of all ranges intersections with the given range.
     *
     * @param r Restriction range.
     */
    virtual void restrictTo(const Range<T>& r) = 0;

    /**
     * @brief Only keep the ranges that fall within the given range.
     *
     * @param r Restriction range.
     */
    virtual void filterWithin(const Range<T>& r) = 0;

    /**
     * @return A string representation of the set of intervals.
     */
    virtual std::string toString() const = 0;

    /**
     * @return True if the set does not contain any range.
     */
    virtual bool isEmpty() const = 0;

    /**
     * @return The number of ranges in the collection.
     */
    virtual size_t size() const = 0;

    /**
     * @return The number of positions in the collection.
     */
    virtual size_t totalLength() const = 0;

    /**
     * @return The ith range in the collection.
     */
    virtual const Range<T>& getRange(size_t i) const = 0;

    /**
     * @brief Clear the collection.
     */
    virtual void clear() = 0;
};

/**
 * @brief A special class used inside RangeCollection.
 */
template<class T> class rangeComp_ {
  public:
    bool operator() (const Range<T>* a, const Range<T>* b) const {
      return ((*a) < (*b));
    }
};

/**
 * @brief This class implements a data structure describing a set of intervales.
 *
 * Intervales can be overlapping, but empty intervales will be ignored/removed.
 */
template<class T> class RangeSet:
  public RangeCollection<T>
{
  public:

  private:
    std::set< Range<T>*, rangeComp_<T> > ranges_;

  public:
    RangeSet(): ranges_() {}

    RangeSet(const RangeSet<T>& set): ranges_()
    {
      for (typename std::set< Range<T>* >::iterator it = set.ranges_.begin(); it != set.ranges_.end(); ++it) {
        ranges_.insert(ranges_.end(), (**it).clone());
      }
    }

    RangeSet& operator=(const RangeSet<T>& set)
    {
      clear_();
      for (typename std::set< Range<T>* >::iterator it = set.ranges_.begin(); it != set.ranges_.end(); ++it) {
        ranges_.insert(ranges_.end(), (**it).clone());
      }
      return *this;
    }

    virtual ~RangeSet() {
      clear_();
    }

  public:
    void addRange(const Range<T>& r) {
      if (!r.isEmpty())
        ranges_.insert(r.clone());
    }

    void restrictTo(const Range<T>& r) {
      typename std::set< Range<T>* >::iterator it = ranges_.begin();
      while (it != ranges_.end()) {
        (**it).sliceWith(r);
        if ((**it).isEmpty()) {
          typename std::set< Range<T>* >::iterator it2 = it;
          delete *it;
          ++it;
          ranges_.erase(it2);
        } else {
          ++it;
        }
      }
    }

    void filterWithin(const Range<T>& r) {
      typename std::set< Range<T>* >::iterator it = ranges_.begin();
      while (it != ranges_.end()) {
        if (r.contains(**it)) {
          ++it;
        } else {
          typename std::set< Range<T>* >::iterator it2 = it;
          delete *it;
          ++it;
          ranges_.erase(it2);
        }
      }
    }
    
    std::string toString() const {
      std::string s = "{ ";
      for (typename std::set< Range<T>* >::const_iterator it = ranges_.begin(); it != ranges_.end(); ++it) {
        s += (**it).toString() + " ";
      }
      s += "}";
      return s;
    }

    bool isEmpty() const { return ranges_.size() == 0; }
    
    size_t size() const { return ranges_.size(); }
    
    /**
     * @return The total length of the range set. Overlapping positions will be counted several times.
     */
    size_t totalLength() const {
      size_t tot = 0;
      for (typename std::set< Range<T>* >::const_iterator it = ranges_.begin(); it != ranges_.end(); ++it) {
        tot += (**it).length();
      }
      return tot;
    }

    const Range<T>& getRange(size_t i) const {
      typename std::set< Range<T>* >::const_iterator it = ranges_.begin();
      for (size_t c = 0; c < i; ++c)
        ++it;
        //it = it++;
      return **it;
    }
    
    const std::set< Range<T>*, rangeComp_<T> >& getSet() const { return ranges_; }

    std::set< Range<T>*, rangeComp_<T> >& getSet() { return ranges_; }

    void clear() {
      clear_();
    }

  private:
    void clear_() {
      for (typename std::set< Range<T>* >::const_iterator it = ranges_.begin(); it != ranges_.end(); ++it) {
        delete *it;
      }
      ranges_.clear();
    }
};

/**
 * @brief This class implements a data structure describing a set of non-overlapping intervales.
 */
template<class T> class MultiRange:
  public RangeCollection<T>
{
  private:
    std::vector<Range<T>*> ranges_;

  public:
    MultiRange(): ranges_() {}

    MultiRange(const MultiRange<T>& mr): ranges_()
    {
      for (size_t i = 0; i < mr.ranges_.size(); ++i)
        ranges_.push_back(mr.ranges_[i]->clone());
    }

    MultiRange& operator=(const MultiRange<T>& mr)
    {
      clear_();
      for (size_t i = 0; i < mr.ranges_.size(); ++i)
        ranges_.push_back(mr.ranges_[i]->clone());
      return *this;
    }

    virtual ~MultiRange() {
      clear_();
    }


  public:
    void addRange(const Range<T>& r) {
      //this is a bit tricky, as many cases can happen. we have to check how many ranges overlap with the new one:
      std::vector<size_t> overlappingPositions;
      for (size_t i = 0; i < ranges_.size(); ++i) {
        if (ranges_[i]->overlap(r))
          overlappingPositions.push_back(i);
      }
      //check if not overlap:
      if (overlappingPositions.size() == 0) {
        //We simply add the new range to the list:
        ranges_.push_back(r.clone());
      } else {
        //We extend the first overlapping element:
        ranges_[overlappingPositions[0]]->expandWith(r);
        //Now we merge all other overlapping ranges, if any:
        for (size_t i = overlappingPositions.size() - 1; i > 0; --i) {
          //Expand first range:
          ranges_[overlappingPositions[0]]->expandWith(*ranges_[overlappingPositions[i]]);
          //Then removes this range:
          delete ranges_[overlappingPositions[i]];
          ranges_.erase(ranges_.begin() + static_cast<ptrdiff_t>(overlappingPositions[i]));
        }
      }
      clean_();
    }

    void restrictTo(const Range<T>& r) {
      for (size_t i = 0; i < ranges_.size(); ++i)
        ranges_[i]->sliceWith(r);
      clean_();
    }

    void filterWithin(const Range<T>& r) {
      typename std::vector< Range<T>* >::iterator it = ranges_.begin();
      while (it != ranges_.end()) {
        if (r.contains(**it)) {
          ++it;
        } else {
          delete *it;
          it = ranges_.erase(it);
        }
      }
    }
 
    /**
     * @return A string representation of the set of intervals.
     */
    std::string toString() const {
      std::string s = "{ ";
      for (size_t i = 0; i < ranges_.size(); ++i)
        s += ranges_[i]->toString() + " ";
      s += "}";
      return s;
    }

    /**
     * @return A vector with all interval bounds.
     */
    std::vector<T> getBounds() const {
      std::vector<T> bounds;
      for (size_t i = 0; i < ranges_.size(); ++i) {
        bounds.push_back(ranges_[i]->begin());
        bounds.push_back(ranges_[i]->end());
      }
      return bounds;
    }

    /**
     * @return True if the set does not contain any range.
     */
    bool isEmpty() const { return ranges_.size() == 0; }

    size_t size() const { return ranges_.size(); }

    size_t totalLength() const {
      size_t tot = 0;
      for (typename std::vector< Range<T>* >::const_iterator it = ranges_.begin(); it != ranges_.end(); ++it) {
        tot += (**it).length();
      }
      return tot;
    }

    const Range<T>& getRange(size_t i) const { return *ranges_[i]; }

    void clear() {
      clear_();
    }

  private:
    void clean_() {
      //Reorder
      rangeComp_<T> comp;
      std::sort(ranges_.begin(), ranges_.end(), comp);
      //Remove empty intervals:
      typename std::vector< Range<T>* >::iterator it = ranges_.begin();
      while (it != ranges_.end()) {
        if ((**it).isEmpty()) {
          delete *it;
          it = ranges_.erase(it);
        } else {
          ++it;
        }
      }
    }
  private:
    void clear_() {
      for (size_t i = 0; i < ranges_.size(); ++i) {
        delete ranges_[i];
      }
      ranges_.clear();
    }

};

} //end of namespace bpp

#endif //_RANGE_H_
