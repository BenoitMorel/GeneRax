#pragma once

extern "C" {
#include <pll.h>
}

enum ParsingErrorType {
  PET_NOERROR = 0,
  PET_FILE_DOES_NOT_EXISTS,
  PET_INVALID_PARENTHESES,
  PET_INVALID_LABEL,
  PET_NOSEMICOLON,
  PET_DOUBLE_BRANCH_LENGTH,
  PET_INVALID_BRANCH_LENGTH,
  PET_POLYTOMY,
  PET_UNROOTED,
  PET_EMPTY_NODE,
  PET_ONLY_ONE_CHILD,
  PET_TOKEN_AFTER_SEMICOLON,
  PET_LAST
};

const char* getParsingErrorName(ParsingErrorType type);
const char* getParsingErrorDiagnostic(ParsingErrorType type);

struct RTreeParsingError {
  ParsingErrorType type;
  unsigned int offset;
};


pll_rtree_t * custom_rtree_parse_newick(const char *s, 
    bool is_file,
    RTreeParsingError *error);

