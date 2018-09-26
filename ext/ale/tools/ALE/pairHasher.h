template <class T>
inline void hash_combine(std::size_t & seed, const T & v)
{
  std::hash<T> hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

namespace std
{
  template<typename S, typename T> struct hash<std::pair<S, T>>
    {
      inline std::size_t operator()(const std::pair<S, T> & v) const
      {
	std::size_t seed = 0;
	::hash_combine(seed, v.first);
	::hash_combine(seed, v.second);
	return seed;
      }
    };
}
