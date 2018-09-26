//
// File: MapTools.h
// Created by: Julien Dutheil
// Created on: Tue May 13 18:16:10 2003
//

/*
Copyright or © or Copr. CNRS, (November 17, 2004)

This software is a computer program whose purpose is to provide utilitary
classes. This file belongs to the Bio++ Project.

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

#ifndef _MAPTOOLS_H_
#define _MAPTOOLS_H_

#include <map>
#include <vector>

namespace bpp
{

/**
 * @brief A fiew tools working on map objects.
 */
class MapTools
{
	public:
	
		/**
		 * @brief Get a vector of all keys in a map.
		 *
		 * @param myMap the map to check.
		 * @return a vector of all keys.
		 */
		template <class Key, class T, class Cmp >
		static std::vector<Key> getKeys(const std::map<Key, T, Cmp> & myMap)
		{
			std::vector<Key> keys;
			for(typename std::map<Key, T>::const_iterator i = myMap.begin(); i != myMap.end(); i++)
      {
				keys.push_back(i->first);
			}
			return keys;
		}
		
		/**
		 * @brief Get a vector of all keys in a map.
		 *
		 * @param myMap the map to check.
		 * @return a vector of all keys.
		 */
		template <class Key, class T >
		static std::vector<Key> getKeys(const std::map<Key, T> & myMap)
		{
			std::vector<Key> keys;
			for(typename std::map<Key, T>::const_iterator i = myMap.begin(); i != myMap.end(); i++)
      {
				keys.push_back(i->first);
			}
			return keys;
		}
		
		/**
		 * @brief Get a vector of all values in a map.
		 *
		 * @param myMap the map to check.
		 * @return a vector of all values.
		 */
		template <class Key, class T, class Cmp >
		static std::vector<T> getValues(const std::map<Key, T, Cmp> & myMap)
		{
			std::vector<T> values;
			for(typename std::map<Key, T>::const_iterator i = myMap.begin(); i != myMap.end(); i++)
      {
				values.push_back(i->second);
			}
			return values;
		}

		/**
		 * @brief Get a vector of all values in a map.
		 *
		 * @param myMap the map to check.
		 * @return a vector of all values.
		 */
		template <class Key, class T >
		static std::vector<T> getValues(const std::map<Key, T> & myMap)
		{
			std::vector<T> values;
			for(typename std::map<Key, T>::const_iterator i = myMap.begin(); i != myMap.end(); i++)
      {
				values.push_back(i->second);
			}
			return values;
		}

};

} //end of namespace bpp.

#endif	//_MAPTOOLS_H_

