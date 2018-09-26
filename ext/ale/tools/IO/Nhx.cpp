// Treerecs – Copyright © by INRIA – All rights reserved – 2017
// Created by Nicolas Comte on 03/09/17.
//

#include "Nhx.h"

#include <ale/Constants.h>

Nhx::Nhx(const bool useTagsAsPptNames) : bpp::Nhx(useTagsAsPptNames) {
  // Create bootstrap property
  //Property bootstrap_property("bootstrap", "B", true, 2);

  // Create species property
  //Property species_property("species", "S", false, 0);

  // Create duplication property
  Property duplication_property(DUPLICATION_STR_FLAG, "D", false, 0);

  // Register properties.
  //registerProperty(bootstrap_property);
  //registerProperty(species_property);
  registerProperty(duplication_property);
}
