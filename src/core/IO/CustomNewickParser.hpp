#pragma once

extern "C" {
#include <pll.h>
}

enum ParsingErrorType {
  PET_NOERROR,
  PET_FILE_DO_NOT_EXISTS,
  PET_INVALID_SYNTAX,
  PET_POLYTOMY
};

struct RTreeParsingError {
  ParsingErrorType type;
};


pll_rtree_t * custom_rtree_parse_newick(const char *s, 
    bool is_file,
    RTreeParsingError *error);

