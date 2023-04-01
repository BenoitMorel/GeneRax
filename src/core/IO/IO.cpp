#include "IO.hpp"
#include <algorithm>

namespace IO {

bool isBlanck(const std::string &s) 
{
  return s.empty() 
    || std::all_of(s.begin(), 
        s.end(), 
        [](char c){return std::isspace(c);});
}

template<typename T, typename P>
T remove_if(T beg, T end, P pred)
{
  T dest = beg;
  for (T itr = beg;itr != end; ++itr)
    if (!pred(*itr))
      *(dest++) = *itr;
  return dest;
}

void removeSpaces(std::string &str) {
  std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
  str.erase(end_pos, str.end());
}



}

