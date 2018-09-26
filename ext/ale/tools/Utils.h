// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 05/12/16.
//

#ifndef PHYLASOLVER_UTILS_H
#define PHYLASOLVER_UTILS_H

//Include std
#include <map>
#include <unordered_map>
#include <vector>
#include <list>
#include <iostream>
#include <iomanip>
#include <cstring>
#include <random>
#include <cassert>

//Include Bpp
#include <Bpp/Phyl/Io/Newick.h>

//Include constants
#include <ale/Constants.h>


// Print in std::ostream stl elements.
template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &pair);

template<class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v);

template<class T>
std::ostream &operator<<(std::ostream &os, const std::list<T> &l);

template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::map<T1, T2> &map);

template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::unordered_map<T1, T2> &map);


// Print in std::ostream bpp elements
std::ostream &operator<<(std::ostream &os, const bpp::PhyloTree &tree);

std::ostream &operator<<(std::ostream &os, const std::shared_ptr<bpp::PhyloNode> &node_ptr);

namespace Utils {
  template<typename T>
  void write_array(std::ostream& os, const T* elements, const std::size_t n){
    os << "{";
    for(std::size_t i = 0; i < n; i++){
      os << elements[i];
      if(i < (n - 1)) os << ", ";
    }
    os << "}";
  }

  template<typename Iterator>
  std::ostream& write(std::ostream &os
      , Iterator it
      , const Iterator& itend
      , const std::string& opening_str = "{"
      , const std::string& elements_separator = ", "
      , const std::string& closing_str = "}")
  {
    os << opening_str;
    if(it != itend) {
      os << *it;
      it++;
      while (it != itend) {
        os << elements_separator << *it;
        it++;
      }
    }
    os << closing_str;

    return os;
  }

  /// Check if a given const char* means "yes" or "true".
  bool isYes(const char* optarg);

  template<typename T>
  T sum(const std::vector<T> &x) {
    T res = T(0.0);
    for (auto &i: x) res += i;
    return res;
    //return std::accumulate(x.begin(), x.end(), 0.0);
  }

  inline double random(const double min = 0.0, const double max = 1.0) {
    std::uniform_real_distribution<double> uniform_real_distribution(min, max);
    return uniform_real_distribution(DEFAULT_RANDOM_GENERATOR);
    // return (max - min) * ( (double)rand() / (double)RAND_MAX ) + min;
    // return bpp::RandomTools::giveRandomNumberBetweenZeroAndEntry(1);
  }

  template<typename T>
  std::vector<std::size_t> getOrder(const std::vector<T> &x, const bool ascending = true) {
    /// Returns a vector of indexes in the order of their corresponding value.
    /// Ex: std::vector<int> x = {3, 1, 2};
    ///     getOrder(x [, true]) will return {1, 2, 0} if ascending = true. 1 is the index of the lowest value in the
    ///     vector and 0 the index of the greatest element.
    ///     getOrder(x, false) will return {0, 2, 1} if ascending = false.

    std::vector<std::size_t> order(x.size());
    std::size_t n(0);
    std::generate(std::begin(order), std::end(order), [&n] { return n++; });
    std::sort(std::begin(order), std::end(order),
              [&x, &ascending](const std::size_t &i, const std::size_t &j) {
                return (ascending) == (x.at(i) < x.at(j));
              }
    );
    return order;
  }

  template<typename Weight>
  std::vector<std::size_t> weightedRandomPick(const std::vector<Weight> &weights, const std::size_t n) {
    /// Random index pick according to the value in weights.
    /// For example, with [0.9, 0.1, 0.0], index 0 will be picked frequently, index 2 never.
    //Init the vector which contains picked indexes
    std::vector<std::size_t> picked;
    picked.reserve(n);

    Weight total = sum(weights);
    auto weights_indexes_in_order = getOrder(weights, true); // Get the descending order of each element of weight

    for (std::size_t i = 0; i < n; i++) {
      //double rand_value = rand() / (double) RAND_MAX;
      double rand_value = RANDOM_DOUBLE_BETWEEN_0_AND_1(DEFAULT_RANDOM_GENERATOR);//random();
      double cumul_prob_value = 0.0;
      for (std::size_t wi: weights_indexes_in_order) {
        cumul_prob_value += ((double) weights.at(wi) / total);
        if (cumul_prob_value > rand_value) {
          picked.push_back(wi);
          break;
        }
      }
    }

    assert(picked.size() == n);

    return picked;
  }

  template<typename Weight>
  std::size_t weightedRandomPick(const std::vector<Weight> &weights) {
    /// Random index pick according to the value in weights.
    /// For example, with [0.9, 0.1, 0.0], index 0 will be picked frequently, index 2 never.
    auto values = weightedRandomPick(weights, 1);
    assert(values.size() == 1);
    return values.at(0);
  }

  template<typename T>
  inline T linearInterpolation(const T& v0, const T& v1, const T& t)
  {
    return (1.0 - t) * v0 + t * v1;
  }

  /*!
   * @brief Get quantiles of a sample of numbers.
   * @tparam Iterator
   * @param begin First element (as an iterator) of the sample to evaluate.
   * @param end Last element (as an iterator) of the sample to evaluate.
   * @param n Type of quantile (2 = median, 3 = tertile, 4 = quartile, ...)
   * @return A vector containing quantiles (as std::vector<double>).
   */
  template<typename Iterator>
  std::vector<double> quantile(
      const Iterator &begin,
      const Iterator &end,
      const std::size_t &n)
  {

    // Change name of template value type for "T".
    typedef typename std::iterator_traits<Iterator>::value_type T;

    std::vector<T> data {begin, end};

    if (data.size() == 0 or data.size() == 1 or n == 1)
    {
      return {data.begin(), data.end()};
    } else if(n < 1) {
      return {};
    }

    std::vector<double> probs(n - 1);

    probs[0] = (1.0/(double)n);

    std::size_t i = 1;
    std::generate(probs.begin() + 1, probs.end(), [&i, p = probs.front()](){ return (p * (++i)); });

    std::sort(data.begin(), data.end());
    std::vector<double> quantiles;
    quantiles.reserve(probs.size());

    for (i = 0; i < probs.size(); ++i)
    {
      double index_li = linearInterpolation<double>(0, data.size() - 1, probs[i]);

      auto left_index = std::max(int64_t(std::floor(index_li)), int64_t(0));
      auto right_index = std::min(int64_t(std::ceil(index_li)), int64_t(data.size() - 1));

      auto left_element = static_cast<double>(data.at(left_index));
      auto right_element = static_cast<double>(data.at(right_index));

      auto quantile = linearInterpolation<double>(left_element, right_element, index_li - left_index);

      quantiles.push_back(quantile);
    }

    return quantiles;
  }

  void printNodeContent(const bpp::PhyloTree &tree, const std::shared_ptr<bpp::PhyloNode> &node,
                        std::ostream &os = std::cout);

  void print_temp(bpp::PhyloTree &tree, std::ostream &os);

  inline bool case_insensitive_char_comp(const unsigned char &a, const unsigned char &b) {
    return std::tolower(a) == std::tolower(b);
  }

  inline bool case_sensitive_char_comp(const unsigned char &a, const unsigned char &b) {
    return a == b;
  }

  bool string_comp(const std::string &a, const std::string &b, const bool case_sensitive = true);

  /// Check if a given string matches with a regex pattern.
  bool strmatch_regex(const std::string &str, const std::string &pattern);

  /// Returns count of non-overlapping occurrences of 'sub' in 'str'.
  unsigned int count(const std::string& str, const std::string& sub);

  std::vector<std::string> splitString(const std::string &str, const char *c, const bool concatenate_delimiters = true);

  template<class T>
  std::size_t getStreamObjectSize(const T &o, const std::ostream &os) {
    /// Returns length in stream of a given element.
    std::ostringstream temp;
    auto old_precision = os.precision();
    temp << std::setprecision(old_precision) << o;
    return temp.str().size();
  }

  /// Return a file name, without path and file extension.
  std::string extractFilename(const std::string& str);

  inline std::string trunc_string_number(const std::string& s, const std::size_t n_dec = 2) {
    auto dot_pos = std::find(s.begin(), s.end(), '.');
    for(std::size_t i = 0; i <= n_dec and dot_pos != s.end(); ++i) dot_pos++;
    return std::string(s.begin(), dot_pos);
  }

  template<typename IterA, typename IterB>
  bool comp_all(IterA a_begin, const IterA &a_end, IterB b_begin, const IterB& b_end) {
    /// Compare two sorted containers.
    if(std::distance(a_begin, a_end) != std::distance(b_begin, b_end))
      return false;

    while(a_begin != a_end){
      if(*a_begin != *b_begin) return false;
      a_begin++;
      b_begin++;
    }
    return true;
  }

  template<typename IterA, typename IterB, typename Comp>
  bool comp_all(IterA a_begin, const IterA &a_end, IterB b_begin, const IterB& b_end, const Comp& comparator) {
    /// Compare two sorted containers using a template comparator.
    if(std::distance(a_begin, a_end) != std::distance(b_begin, b_end))
      return false;

    while(a_begin != a_end){
      if(not comparator(*a_begin, *b_begin)) return false;
      a_begin++;
      b_begin++;
    }
    return true;
  }

  /// Returns a list of value indexes which have a match with a template condition in a container.
  /// For example:
  ///
  ///   std::vector<int> my_vect {1, 42, 3, 300, 28};
  ///   std::cout << "Values > 30 in my_vect:" ;
  ///   auto indexes = matchesValuesIndexes(my_vect.begin(), my_vect.end(), [](int x){ return x > 30; });
  ///   for(auto i: indexes) std::cout << " " << my_vect.at(i);
  ///   std::cout << "." << std::endl;
  ///
  /// \tparam Iterator
  /// \tparam Condition
  /// \param begin
  /// \param end
  /// \param condition
  /// \return std::list of indexes.
  template<typename Iterator, typename Condition>
  std::list<std::size_t> matchesValueIndexes(
      Iterator begin, const Iterator& end, const Condition& condition) {
    std::list<std::size_t> res;
    std::size_t i = 0;
    while(begin != end){
      if(condition(*begin)) {
        res.push_back(i);
      }
      i++;
      begin++;
    }
    return res;
  }

  /// Check if a container has an element which match with a condition
  template<typename Iteration, typename Lambda>
  bool contains(const Iteration& begin, const Iteration& end, const Lambda& lambda){
    auto current = begin;
    while(current != end){
      if(lambda(*current)) return true;
      current++;
    }
    return false;
  }


  /// Trim at left.
  inline void ltrim_str(std::string &str) {
    str.erase(str.begin(), std::find_if(str.begin(), str.end(), [](int ch) {
      return !std::isspace(ch);
    }));
  }

  /// Trim at right.
  inline void rtrim_str(std::string &str) {
    str.erase(std::find_if(str.rbegin(), str.rend(), [](int ch) {
      return !std::isspace(ch);
    }).base(), str.end());
  }

  /// Trim at left and right.
  inline void trim_str(std::string &str) {
    Utils::ltrim_str(str);
    Utils::rtrim_str(str);
  }

  /// Trim at right.
  inline std::string ltrim(std::string str) {
    Utils::ltrim_str(str);
    return str;
  }

  /// Trim at left and right.
  inline std::string rtrim(std::string str) {
    Utils::rtrim_str(str);
    return str;
  }

  /// Trim at left and right.
  inline std::string trim(std::string str) {
    Utils::trim_str(str);
    return str;
  }

  /// Returns a string with a specific substring replaced by an other
  /// \param str
  /// \param old_substr Substring to replace, can be a regex pattern.
  /// \param new_substr New element to put in place of the previous substring/pattern.
  /// \return std::string with substrings replaced (all occurrences).
  std::string replace(const std::string& str
                      , const std::string& old_substr
                      , const std::string& new_substr);

  /// Returns a string with a specific char replaced by an other
  /// \param str
  /// \param old_char character to replace.
  /// \param new_char New element to put in place of the previous character.
  /// \return std::string with substrings replaced (all occurrences).
  std::string replace(const std::string& str
                      , const char old_char
                      , const char new_char);

  /// Check if a string has a match with a number.
  bool stringIsNumber(const std::string &s);

  /// Compare two doubles according to a specific precision, defined by delta (equivalent to operator==).
  inline bool double_equivalence(const double a, const double b, const double delta = DEFAULT_DOUBLE_EQUIVALENCE_PRECISION)
  { return fabs(a - b) < delta; }

  /// Compare two doubles (operator <=).
  inline bool double_equal_or_inferior(const double a, const double b) { return double_equivalence(a, b) or (a < b); }

  /// Compare two doubles (operator >=).
  inline bool double_equal_or_superior(const double a, const double b) { return double_equivalence(a, b) or (a > b); }

  /// Factorial(n).
  template<typename T>
  T factorial(T n){
    T res = 1.0;

    for (T i = 1.0 ; i <= n ; ++i) { // warning: double precision cannot be fine with equality
      res *= i;
    }

    return res;
  }

  /// Double factorial function.
  template<typename T>
  T double_factorial(T n) {
    T res = 1.0;

    for(T i = n ; i >= 1.0 ; i-= 2.0){ // warning: double precision cannot be fine with equality
      res *=i;
    }

    return res;
  }
}

template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::pair<T1, T2> &pair) {
  /// Print a std::pair.
  os << "(" << pair.first << ", " << pair.second << ")";
  return os;
}

/*!
 * @brief Prints each element of a given std::vector.
 * @param os
 * @param v The vector that contains each element to print
 * @return std::ostream&
 */
template<class T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &v) {
  ///Streams all elements of a vector (std).
  return Utils::write(os, v.begin(), v.end(), "[", ", ", "]");
}

template<class T>
std::ostream &operator<<(std::ostream &os, const std::list<T> &l) {
  /// Streams all elements of a list (std).
  return Utils::write(os, l.begin(), l.end(), "[", ", ", "]");
}

template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::map<T1, T2> &map) {
  /// Streams all pairs of a map (std).
  return Utils::write(os, map.begin(), map.end(), "{", ", ", "}");
}

template<class T1, class T2>
std::ostream &operator<<(std::ostream &os, const std::unordered_map<T1, T2> &map) {
  /// Streams all pairs of an unordered_map (std).
  return Utils::write(os, map.begin(), map.end(), "{", ", ", "}");
}

std::ostream &operator<<(std::ostream &os, const bpp::PhyloTree &tree);

std::ostream &operator<<(std::ostream &os, const std::shared_ptr<bpp::PhyloNode> &node_ptr);

#endif //PHYLASOLVER_UTILS_H
