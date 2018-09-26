// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 03/09/17.
//

#ifndef TREERECS_NHX_H
#define TREERECS_NHX_H

#include <Bpp/Phyl/Io/Nhx.h>

/*!
 * @class Nhx
 * @brief Nhx provides functions to load a bpp::PhyloTree or write it into a Nhx file (inheritance from bpp::Nhx).
 */
class Nhx: public bpp::Nhx {
public:
  explicit Nhx(const bool useTagsAsPptNames = false);
};


#endif //TREERECS_NHX_H
