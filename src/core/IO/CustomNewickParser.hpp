#pragma once

extern "C" {
#include <pll.h>
}

enum ParsingErrorType {
  PET_NOERROR,
  PET_FILE_DO_NOT_EXISTS,
  PET_INVALID_PARENTHESIS,
  PET_INVALID_LABEL,
  PET_INVALID_SYNTAX,
  PET_NOSEMICOLON,
  PET_DOUBLE_BRANCH_LENGTH,
  PET_INVALID_BRANCH_LENGTH,
  PET_POLYTOMY,
  PET_EMPTY_NODE,
  PET_ONLY_ONE_CHILD,
  PET_TOKEN_AFTER_SEMICOLON
};

struct RTreeParsingError {
  ParsingErrorType type;
};


pll_rtree_t * custom_rtree_parse_newick(const char *s, 
    bool is_file,
    RTreeParsingError *error);

