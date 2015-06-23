#ifndef UTILS_H
#define UTILS_H

#include <algorithm>
#include <tr1/functional>

template<class RandIt, class Compare>
  bool next_k_permutation(RandIt first, RandIt mid, RandIt last, Compare comp)
{
  std::sort(mid, last, std::tr1::bind(comp, std::tr1::placeholders::_2
				      , std::tr1::placeholders::_1));
  return std::next_permutation(first, last, comp);
}

template<class RandIt>
bool next_k_permutation(RandIt first, RandIt mid, RandIt last)
{
  typedef typename std::iterator_traits<RandIt>::value_type value_type;
  
  std::sort(mid, last, std::greater< value_type >());
  return std::next_permutation(first, last);
}

template<class BiDiIt, class Compare>
  bool next_combination(BiDiIt first, BiDiIt mid, BiDiIt last, Compare comp)
{
  bool result;
  do
    {
      result = next_k_permutation(first, mid, last, comp);
    } while (std::adjacent_find( first, mid,
				 std::tr1::bind(comp, std::tr1::placeholders::_2
						, std::tr1::placeholders::_1) )
	     != mid );
  return result;
}

#endif
