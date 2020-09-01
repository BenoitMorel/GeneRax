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

/**
 *  Get a short name describing the parsing error
 */
const char* getParsingErrorName(ParsingErrorType type);

/**
 * Get a short description of the parsing error
 */
const char* getParsingErrorDiagnostic(ParsingErrorType type);

/*
 *  Parsing error information
 */
struct ParsingError {
  // Error type. See also getParsingErrorName and
  // getParsingErrorDiagnostic
  ParsingErrorType type;

  // Offset in the newick string at which the error
  // was detected
  unsigned int offset;
};

/**
 *  Parse a rooted tree from a newick string or file
 *  If an error occures, returns NULL and fills the 
 *  error object.
 */
pll_rtree_t * custom_rtree_parse_newick(const char *s, 
    bool is_file,
    ParsingError *error);


