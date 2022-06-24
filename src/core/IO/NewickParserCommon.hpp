#pragma once

#include <corax/corax.h>
#include <IO/LibpllParsers.hpp>

extern int is_separator[];
extern int to_trim[];


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
const char* get_parsing_error_name(ParsingErrorType type);

/**
 * Get a short description of the parsing error
 */
const char* get_parsing_error_diagnostic(ParsingErrorType type);

/*
 *  Parsing error information
 */
struct ParsingError {
  // Error type. See also get_parsing_error_name and
  // get_parsing_error_diagnostic
  ParsingErrorType type;

  // Offset in the newick string at which the error
  // was detected
  unsigned int offset;
};


int is_numeric(char *s, double *d);

/**
 *  Increment buffer until next character is not
 *  a space or tab
 */
void skip_spaces(char **buffer);


/**
 *  Return the size of the next token to read 
 *  (number of characters before the next separator
 *  after the current character)
 */
unsigned int get_token_size(char *buffer);

/**
 *  Reads and return the content of a file.
 *  Deallocation has to be done by the caller with free
 */
char *get_file_content(const char *filename);

enum TokenType {
  TT_LEFT_PAR,
  TT_RIGHT_PAR,
  TT_COMMA,
  TT_COLON,
  TT_SEMICOLON,
  TT_DOUBLE,
  TT_STRING
};

struct Token {
  // type of token (string, float, special character etc.)
  TokenType type;
  // value of the token (only relevant for string and float tokens)
  char *str;
  // size of str
  unsigned int str_size;
  // floating value of the token (only relevant for float tokens)
  double numeric_value;
};

/**
 *  Fill token with the next token in buffer,
 *  and increment the buffer pointer to the next pointer
 *  The token should then be destroyed with destroy_token
 */
unsigned int read_token(char **buffer, 
    Token *token);

/**
 *  Deallocate the elements of token created
 *  with read_token
 */
void destroy_token(Token *token);




