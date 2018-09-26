//
// File GraphObserver.h
// Created by: Thomas Bigot
// Last modification : vendredi 4 novembre 2016, à 10h 21
//

/*
   Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

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

#ifndef _GRAPHOBSERVER_HPP_
#define _GRAPHOBSERVER_HPP_

#include "Graph.h"

#include "../Exceptions.h"
#include "../Clonable.h"

#include <vector>
#include <map>
#include <iostream>
#include <ostream>


namespace bpp
{
/**
 * @brief Defines a Graph Observer. It is a template which follows
 * (subscribed to) a Graph.
 * The graph and the graph observer communicate to keep them up-to-date
 * each other. The observer is also an actor, since it can change
 * the structure of the observed Graph.
 *
 * @author Thomas Bigot
 */

// interface

class GraphObserver :
  public virtual Clonable
{
public:
  /** @name Function called by the subjectGraph
   *  Methods called by the subject graph to make this observer so fit the subject graph
   */
  // /@{

  /**
   * Delete unused object edges, since they have been deleted in the graph
   * @param edgesToDelete a vector of Edges to delete
   */
  virtual void deletedEdgesUpdate(const std::vector< unsigned int >& edgesToDelete) = 0;

  /**
   * Delete unused object nodes, since they have been deleted in the graph
   * @param nodesToDelete a vector of N to delete
   */
  virtual void deletedNodesUpdate(const std::vector< unsigned int >& nodesToDelete) = 0;

  // /@}
};
}

#else

namespace bpp
{class GraphObserver; }

#endif
