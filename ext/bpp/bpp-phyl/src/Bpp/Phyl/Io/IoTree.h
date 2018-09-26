//
// File: IOTree.h
// Created by: Julien Dutheil
// Created on: Thu Oct 23 15:19:06 2003
//

/*
  Copyright or © or Copr. CNRS, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

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

#ifndef _IOTREE_H_
#define _IOTREE_H_

#include "../Tree/PhyloTree.h"

// From the STL:
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>

#include <Bpp/Exceptions.h>
#include <Bpp/Io/IoFormat.h>

namespace bpp
{

  /**
   * @brief General interface for tree I/O.
   */
  class IOTree:
    public virtual IOFormat
  {
  public:
    IOTree() {}
    virtual ~IOTree() {}

  public:
    virtual const std::string getDataType() const { return "Tree"; }
  };

/**
 * @brief General interface for tree readers.
 */
  class ITree:
    public virtual IOTree
  {
  public:
    ITree() {}
    virtual ~ITree() {}

  public:
    /**
     * @brief Read a tree from a file.
     *
     * @param path The file path.
     * @return A new tree object.
     * @throw Exception If an error occured.
     */

    virtual PhyloTree* readP(const std::string& path) const = 0;
    
    /**
     * @brief Read a tree from a stream.
     *
     * @param in The input stream.
     * @return A new tree object.
     * @throw Exception If an error occured.
     */

    virtual PhyloTree* readP(std::istream& in) const = 0;
  };

/**
 * @brief General interface for tree writers.
 */
  class OTree:
    public virtual IOTree
  {
  public:
    OTree() {}
    virtual ~OTree() {}

  public:
    /**
     * @brief Write a tree to a file.
     *
     * @param tree A tree object.
     * @param path The file path.
     * @param overwrite Tell if existing file must be overwritten.
     * Otherwise append to the file.
     * @throw Exception If an error occured.
     */
    virtual void write(const PhyloTree& tree, const std::string& path, bool overwrite) const = 0;
/**
     * @brief Write a tree to a stream.
     *
     * @param tree A tree object.
     * @param out The output stream.
     * @throw Exception If an error occured.
     */
    virtual void write(const PhyloTree& tree, std::ostream& out) const =  0;

  };

/**
 * @brief Partial implementation of the ITree interface.
 */
  class AbstractITree:
    public virtual ITree
  {

    /*
     * @brief Basic element for branch description.
     *
     */
    
  protected:
    struct Element
    {
    public:
      std::string content;
      std::string length;
      std::string annotation;
      bool isLeaf;

    public:
      Element() : content(),
                  length(),
                  annotation(),
                  isLeaf(false) {}
    };

  public:
    AbstractITree() {}
    virtual ~AbstractITree() {}

  public:
    virtual PhyloTree* readP(std::istream& in) const = 0;
    virtual PhyloTree* readP(const std::string& path) const
    {
      std::ifstream input(path.c_str(), std::ios::in);
      PhyloTree* tree = readP(input);
      input.close();
      return tree;
    }

    virtual Element getElement(const std::string& elt)
    {
      return Element();
    }
    
  };

/**
 * @brief Partial implementation of the OTree interface.
 */
  class AbstractOTree:
    public virtual OTree
  {
  public:
    AbstractOTree() {}
    virtual ~AbstractOTree() {}

  public:
    void write(const PhyloTree& tree, std::ostream& out) const = 0;
    virtual void write(const PhyloTree& tree, const std::string& path, bool overwrite) const
    {
      try {
        // Open file in specified mode

        std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out|std::ios::app));
        write(tree, output);
        output.close();
      }
      catch (IOException e)
      {
        std::stringstream ss ;
        ss << e.what() <<"\nProblem writing tree to file "<< path <<"\n Is the file path correct and do \
you have the proper authorizations? ";
        throw (IOException ( ss.str() ) );
      }

    }
  };




/**
 * @brief General interface for multiple trees readers.
 */
  class IMultiTree:
    public virtual IOTree
  {
  public:
    IMultiTree() {}
    virtual ~IMultiTree() {}

  public:
    /**
     * @brief Read trees from a file.
     *
     * @param path The file path.
     * @param trees The output trees container.
     * @throw Exception If an error occured.
     */
    virtual void read(const std::string& path, std::vector<PhyloTree*>& trees) const = 0;
    /**
     * @brief Read trees from a stream.
     *
     * @param in The input stream.
     * @param trees The output trees container.
     * @throw Exception If an error occured.
     */
    virtual void read(std::istream& in, std::vector<PhyloTree*>& trees) const = 0;
  };

/**
 * @brief General interface for tree writers.
 */
  class OMultiTree:
    public virtual IOTree
  {
  public:
    OMultiTree() {}
    virtual ~OMultiTree() {}

  public:
    /**
     * @brief Write trees to a file.
     *
     * @param trees A vector of tree objects.
     * @param path The file path.
     * @param overwrite Tell if existing file must be overwritten.
     * Otherwise append to the file.
     * @throw Exception If an error occured.
     */

    virtual void write(const std::vector<const PhyloTree*>& trees, const std::string& path, bool overwrite) const = 0;

    /**
     * @brief Write trees to a stream.
     *
     * @param trees A vector of tree objects.
     * @param out The output stream.
     * @throw Exception If an error occured.
     */

    virtual void write(const std::vector<const PhyloTree*>& trees, std::ostream& out) const = 0;
  };

  /**
   * @brief Partial implementation of the IMultiTree interface.
   */
  class AbstractIMultiTree:
    public virtual IMultiTree
  {
  public:
    AbstractIMultiTree() {}
    virtual ~AbstractIMultiTree() {}

  public:
    virtual void read(std::istream& in, std::vector<PhyloTree*>& trees) const = 0;
    virtual void read(const std::string& path, std::vector<PhyloTree*>& trees) const
    {
      std::ifstream input(path.c_str(), std::ios::in);
      read(input, trees);
      input.close();
    }

  };

/**
 * @brief Partial implementation of the OTree interface.
 */
  class AbstractOMultiTree:
    public virtual OMultiTree
  {
  public:
    AbstractOMultiTree() {}
    virtual ~AbstractOMultiTree() {}

  public:
    void write(const std::vector<const PhyloTree*>& trees, std::ostream& out) const = 0;
    virtual void write(const std::vector<const PhyloTree*>& trees, const std::string& path, bool overwrite) const

    {
      // Open file in specified mode
      std::ofstream output(path.c_str(), overwrite ? (std::ios::out) : (std::ios::out|std::ios::app));
      write(trees, output);
      output.close();
    }
  };

} //end of namespace bpp.

#endif  //_IOTREE_H_

