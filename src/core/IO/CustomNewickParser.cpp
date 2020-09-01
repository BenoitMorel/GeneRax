#include "CustomNewickParser.hpp"
#include <ctype.h>
#include <stdlib.h>
#include <iostream>

int is_separator[256] = {1,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
int to_trim[256] = {0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


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

const char *error_type_diagnostic[PET_LAST] = {
  "Tree was successfully parsed",
  "File does not exist",
  "The newick string contains too few left or right parentheses, or a semicolon was inserted too early.",
  "A node label is invalid or was set twice. Please also check the presence of invalid characters, like ()[]'\";,: newlines spaces and tabs",
  "The newick string should end with a semicolon",
  "The branch length of a node was set twice",
  "The branch length of a node could not be read as a floating value",
  "The tree contains polytomies (nodes with strictly more than two children",
  "The tree is unrooted: its top node has strictly more than two children",
  "A node has no label nor children",
  "A node has only one child",
  "Some text follows the ending semicolon"
};


const char* getParsingErrorName(ParsingErrorType type)
{
  return error_type_name[type];
}

const char* getParsingErrorDiagnostic(ParsingErrorType type)
{
  return error_type_diagnostic[type];
}

char *getFileContent(const char *filename)
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

struct RTreeParser {
  char* input;
  char* input_current;
  unsigned int input_size;
  pll_rnode_t **nodes;
  unsigned int nodes_capacity;
  unsigned int nodes_number;
  pll_rnode_t *current_node;
  pll_rnode_t *parent_node;
  bool is_file;
  RTreeParsingError *error;
};

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
  TokenType type;
  char *str;
  unsigned int str_size;
  double numeric_value;
};

int has_errored(struct RTreeParser *p)
{
  return p->error->type != PET_NOERROR;
}

void set_error(struct RTreeParser *p, ParsingErrorType type)
{
  if (!has_errored(p)) {
    p->error->type = type;
    p->error->offset = p->input_current - p->input - 1;
  }
}

void destroy_rtree_parser(struct RTreeParser *p)
{
  for (unsigned int i = 0; i < p->nodes_number; ++i) {
    if (p->nodes[i]) {
      free(p->nodes[i]->label);
    }
    free(p->nodes[i]);
  }
  free(p->nodes);
  if (p->is_file) {
    free(p->input);
  }
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

unsigned int get_string_size(RTreeParser *p)
{
  char *curr = p->input_current;
  while (!is_separator[(int)*curr]) {
    curr++;
  }
  return curr - p->input_current; 
}

unsigned int read_token(RTreeParser *p, 
    Token *token)
{
  char **buffer = &p->input_current;
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
    unsigned int string_size = get_string_size(p);
    token->str = (char *)malloc(string_size + 1);
    memcpy(token->str, *buffer, string_size);
    token->str[string_size] = '\0';
    token->type = is_numeric(token->str, &token->numeric_value)
      ? TT_DOUBLE : TT_STRING;
    *buffer += string_size;
  }
  return 1;
}

void increase_nodes_capacity(RTreeParser *p)
{
  p->nodes_capacity *= 4;  
  pll_rnode_t **new_buffer = (pll_rnode_t **)calloc(
      p->nodes_capacity, sizeof(pll_rnode_t *));
  memcpy(new_buffer, p->nodes, (p->nodes_number) * sizeof(pll_rnode_t *));
  free(p->nodes);
  p->nodes = new_buffer;
}

pll_rnode_t *rtree_parse_add_node(RTreeParser *p)
{
  if (p->nodes_number >= p->nodes_capacity) {
    increase_nodes_capacity(p);
  }
  pll_rnode_t *node = (pll_rnode_t *)malloc(sizeof(pll_rnode_t));
  p->nodes[p->nodes_number] = node;
  node->node_index = p->nodes_number;
  node->label = NULL;
  node->length = INFINITY;
  node->left = NULL;
  node->right = NULL;
  node->data = NULL;
  node->parent = p->parent_node;
  if (node->parent != NULL) {
    pll_rnode_t *parent_node = node->parent;
    if (!parent_node->left) {
      parent_node->left = node;
    } else if (!parent_node->right) {
      parent_node->right = node;
    } else {
      if (parent_node->parent == NULL) {
        set_error(p, PET_UNROOTED);
      } else {
        set_error(p, PET_POLYTOMY);
      }
      free(node);
      p->nodes[p->nodes_number] = NULL;
      return NULL;
    }
  }
  p->nodes_number++;
  return node; 
}

void rtree_add_node_down(RTreeParser *p) 
{
  p->parent_node = p->current_node;
  p->current_node = rtree_parse_add_node(p);
}

void rtree_go_up(RTreeParser *p)
{
  if (!p->current_node->label && !p->current_node->left) {
    set_error(p, PET_EMPTY_NODE);
  }
  if (p->current_node->left && !p->current_node->right) {
    set_error(p, PET_ONLY_ONE_CHILD);
  }
  p->current_node = p->parent_node;
  p->parent_node = p->parent_node->parent;
}

      
void rtree_add_node_neighbor(RTreeParser *p)
{
  p->current_node = rtree_parse_add_node(p);
}

int is_branch_length_set(pll_rnode_t *node)
{
  return node && node->length != INFINITY;
}

void rtree_add_label(RTreeParser *p, Token *token)
{
  if (p->current_node->label) {
    set_error(p, PET_INVALID_LABEL);
    return;
  }
  if (is_branch_length_set(p->current_node)) {
    set_error(p, PET_INVALID_BRANCH_LENGTH);
  }
  p->current_node->label = token->str;
  token->str = NULL;
}

void destroy_token(Token *token)
{
  free(token->str);
  token->str = NULL;
}


void terminate_node_creation(pll_rnode_t *node)
{
  if (!node) {
    return;
  }
  if (!is_branch_length_set(node)) {
    node->length = 0.0;
  }
}

void parse(RTreeParser *p)
{
  Token token;
  Token tokenBL;
  token.str = NULL;
  tokenBL.str = NULL;
  bool end = false;
  rtree_add_node_down(p); // add root
   
  while (!end && read_token(p, &token)) 
  {
    switch(token.type) {
    case TT_SEMICOLON:
      end = true;
      break;
    case TT_COLON:
      read_token(p, &tokenBL);
      if (tokenBL.type != TT_DOUBLE) {
        set_error(p, PET_INVALID_BRANCH_LENGTH);
      }
      if (is_branch_length_set(p->current_node)) {
        set_error(p, PET_DOUBLE_BRANCH_LENGTH);
      }
      p->current_node->length = tokenBL.numeric_value;
      destroy_token(&tokenBL);
      break;
    case TT_LEFT_PAR:
      // go down in the tree
      rtree_add_node_down(p);
      break;
    case TT_RIGHT_PAR:
      if (!p->parent_node) {
        // at this point, we have more right than
        // left parenthesis
        set_error(p, PET_INVALID_PARENTHESES);
      } else {
        terminate_node_creation(p->current_node);
        rtree_go_up(p); 
      }
      break;
    case TT_COMMA:
      terminate_node_creation(p->current_node);
      rtree_add_node_neighbor(p);
      break;
    case TT_STRING:
    case TT_DOUBLE:
      rtree_add_label(p, &token); 
      break;
    }
    destroy_token(&token);
    end |= has_errored(p);
  }
  if (p->parent_node) {
    // we have more left than right parenthesis
    set_error(p, PET_INVALID_PARENTHESES);
  }
  if (!end) {
    set_error(p, PET_NOSEMICOLON);
  }
  if (read_token(p, &token)) {
    set_error(p, PET_TOKEN_AFTER_SEMICOLON);
    destroy_token(&token);
  }
  if (p->error->type == PET_NOERROR) {
    terminate_node_creation(p->current_node);
  }
}

pll_rtree_t *build_rtree(RTreeParser *p)
{
  if (has_errored(p)) {
    return NULL;
  }
  pll_rtree_t * tree = (pll_rtree_t *)malloc(sizeof(pll_rtree_t));
  tree->nodes = (pll_rnode_t **)malloc(
      (p->nodes_number)*sizeof(pll_rnode_t *));
  unsigned int tips_number = p->nodes_number / 2 + 1;
  unsigned int tips_index = 0;
  unsigned int internal_index = tips_number;
  // start at one, because the first node is the root
  // and should be placed at the end
  for (unsigned int i = 1; i < p->nodes_number; ++i) {
    pll_rnode_t *node = p->nodes[i];
    p->nodes[i] = NULL; // avoid double free when destroying p
    if (!node->left) {
      node->node_index 
        = node->clv_index 
        = node->pmatrix_index 
        = tips_index++;
      node->scaler_index = PLL_SCALE_BUFFER_NONE;
    } else {
      node->scaler_index = internal_index - tips_number;
      node->node_index 
        = node->clv_index 
        = node->pmatrix_index 
        = internal_index++;
    }
    tree->nodes[node->node_index] = node;
  }
  tree->root = p->nodes[0];
  p->nodes[0] = NULL; // avoid double free when destroying p
  tree->root->scaler_index = internal_index - tips_number;
  tree->root->node_index 
        = tree->root->clv_index 
        = internal_index++;
  tree->root->pmatrix_index = 0;
  tree->nodes[tree->root->node_index] = tree->root;
  assert(internal_index == p->nodes_number);
  assert(tips_index == tips_number);
  tree->tip_count = tips_number;
  tree->inner_count = p->nodes_number - tips_number;
  tree->edge_count = p->nodes_number - 1;
  return tree;
}

pll_rtree_t * custom_rtree_parse_newick(const char *input,
    bool is_file,
    RTreeParsingError *error)
{
  RTreeParser p;
  p.error = error;
  error->type = PET_NOERROR;
  // todo: use const correctly!!
  p.is_file = is_file;
  if (is_file) {
    p.input = getFileContent(input); 
    if (p.input == NULL) {
      set_error(&p, PET_FILE_DOES_NOT_EXISTS);
      return NULL;
    }
  } else {
    p.input = (char*)input;
  }
  p.input_current = (char*)p.input;
  p.nodes_capacity = 1000; 
  p.nodes_number = 0;
  p.nodes = (pll_rnode_t **)
    calloc(p.nodes_capacity, sizeof(pll_rnode_t *));
  assert(p.nodes);
  p.current_node = NULL;
  p.parent_node = NULL;
  parse(&p); 
  pll_rtree_t *rtree = build_rtree(&p);
  destroy_rtree_parser(&p);
  return rtree;
}


