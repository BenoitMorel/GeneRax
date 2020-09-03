#include "NewickParserCommon.hpp"
#include <ctype.h>
#include <stdlib.h>
#include <iostream>

// a separator is a character that ends a label
int is_separator[256] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,0,0,0,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};

// characters to trim (outside a label)
int to_trim[256] = {0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


// error type to error name
const char *error_type_name[PET_LAST] = {
  "NOERROR",
  "FILE_DOES_NOT_EXISTS",
  "INVALID_PARENTHESES",
  "INVALID_LABEL",
  "NO_SEMICOLON",
  "DOUBLE_BRANCH_LENGTH",
  "INVALID_BRANCH_LENGTH",
  "POLYTOMY",
  "UNROOTED",
  "EMPTY_NODE",
  "ONLY_ONE_CHILD",
  "TOKEN_AFTER_SEMICOLON"
};


// error type to error diagnostic
const char *error_type_diagnostic[PET_LAST] = {
  "Tree was successfully parsed",
  "File does not exist",
  "The newick string contains too few left or right parentheses, or a semicolon was inserted too early.",
  "A node label is invalid or was set twice. Please also check the presence of invalid characters, like ()[]'\";,: newlines spaces, tabs and unprintable characters like \x01",
  "The newick string should end with a semicolon",
  "The branch length of a node was set twice",
  "The branch length of a node could not be read as a floating value",
  "The tree contains polytomies (nodes with strictly more than two children. This error might also be raised when a comma is missing",
  "The tree is unrooted: its top node has strictly more than two children",
  "A node has no label nor children",
  "A node has only one child",
  "Some text follows the ending semicolon"
};


const char* get_parsing_error_name(ParsingErrorType type)
{
  return error_type_name[type];
}

const char* get_parsing_error_diagnostic(ParsingErrorType type)
{
  return error_type_diagnostic[type];
}

int is_numeric(char *s, double *d)
{
  if (s == NULL || *s == '\0' || isspace(*s))
    return 0;
  char * p;
  *d = strtod (s, &p);
  return *p == '\0';
}


void skip_spaces(char **buffer) 
{
  while (to_trim[(int)**buffer]) {
    (*buffer)++;
  }
}

unsigned int get_token_size(char *buffer)
{
  char *curr = buffer;
  while (!is_separator[(int)*curr]) {
    curr++;
  }
  return curr - buffer; 
}

char *get_file_content(const char *filename)
{
  char * buffer = 0;
  long length;
  FILE * f = fopen (filename, "rb");

  if (f)
  {
    fseek (f, 0, SEEK_END);
    length = ftell (f);
    fseek (f, 0, SEEK_SET);
    buffer = (char *)malloc(length + 1);
    if (buffer) {
      fread (buffer, 1, length, f);
    }
    buffer[length] = '\0';
    fclose (f);
  }
  return buffer;
}

unsigned int read_token(char **buffer, 
    Token *token)
{
  skip_spaces(buffer);
  if (**buffer == '\0') {
    return 0;
  } else if (**buffer == '(') {
    token->type = TT_LEFT_PAR;
    (*buffer)++;
  } else if (**buffer == ')') {
    token->type = TT_RIGHT_PAR;
    (*buffer)++;
  } else if (**buffer == ',') {
    token->type = TT_COMMA;
    (*buffer)++;
  } else if (**buffer == ':') {
    token->type = TT_COLON;
    (*buffer)++;
  } else if (**buffer == ';') {
    token->type = TT_SEMICOLON;
    (*buffer)++;
  } else {
    unsigned int string_size = get_token_size(*buffer);
    token->str = (char *)malloc(string_size + 1);
    memcpy(token->str, *buffer, string_size);
    token->str[string_size] = '\0';
    token->type = is_numeric(token->str, &token->numeric_value)
      ? TT_DOUBLE : TT_STRING;
    *buffer += string_size;
  }
  return 1;
}

void destroy_token(Token *token)
{
  free(token->str);
  token->str = NULL;
}

