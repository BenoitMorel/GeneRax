#include "CustomNewickParser.hpp"
#include <ctype.h>
#include <stdlib.h>
#include <iostream>
static const unsigned int INVALID_INDEX = (unsigned int)-1;

int is_separator[256] = {1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};


#define MAX_BUFFER_SIZE 8824

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
  unsigned int current_node_index;
  unsigned int parent_node_index;
  unsigned int depth;
  bool is_file;
  RTreeParsingError *error;
  char file_buffer[MAX_BUFFER_SIZE];
  char previous_file_buffer[MAX_BUFFER_SIZE];
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

void destroy_rtree_parser(struct RTreeParser *p)
{
  free(p->nodes);
}

int is_numeric(char *s, double *d)
{
  if (s == NULL || *s == '\0' || isspace(*s))
    return 0;
  char * p;
  *d = strtod (s, &p);
  return *p == '\0';
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

  unsigned int capacity = p->nodes_capacity;
  p->nodes_capacity *= 4;  
  pll_rnode_t **new_buffer = (pll_rnode_t **)calloc(
      p->nodes_capacity, sizeof(pll_rnode_t *));
  memcpy(new_buffer, p->nodes, capacity * sizeof(pll_rnode_t *));
  free(p->nodes);
  p->nodes = new_buffer;
}

pll_rnode_t *rtree_parse_add_node(RTreeParser *p)
{
  p->current_node_index = p->nodes_number;
  if (p->current_node_index >= p->nodes_capacity) {
    increase_nodes_capacity(p);
  }
  pll_rnode_t *node = (pll_rnode_t *)malloc(sizeof(pll_rnode_t));
  p->nodes[p->current_node_index] = node;
  node->label = NULL;
  node->length = 0.0;
  node->node_index = p->current_node_index;
  node->left = NULL;
  node->right = NULL;
  node->data = NULL;
  if (INVALID_INDEX != p->parent_node_index) {
    node->parent = p->nodes[p->parent_node_index];
  } else {
    node->parent = NULL;
  }
  if (node->parent != NULL) {
    pll_rnode_t *parent_node = node->parent;
    if (!parent_node->left) {
      parent_node->left = node;
    } else if (!parent_node->right) {
      parent_node->right = node;
    } else {
      assert(false); // polytomy
    }
  }
  p->nodes_number++;
  return node; 
}

void rtree_add_node_down(RTreeParser *p) 
{
  p->depth++;
  p->parent_node_index = p->current_node_index;
  pll_rnode_t *current_node = rtree_parse_add_node(p);
  p->current_node_index = current_node->node_index;
}

void rtree_go_up(RTreeParser *p)
{
  p->depth--;
  if (p->parent_node_index != INVALID_INDEX) {
    auto parent = p->nodes[p->parent_node_index];
  }
  pll_rnode_t *new_current_node = p->nodes[p->parent_node_index];
  p->current_node_index = new_current_node->node_index;
  
  p->parent_node_index = new_current_node->parent ? new_current_node->parent->node_index : INVALID_INDEX;
}

      
void rtree_add_node_neighbor(RTreeParser *p)
{
  pll_rnode_t *current_node = rtree_parse_add_node(p);
  p->current_node_index = current_node->node_index;
}

void rtree_add_label(RTreeParser *p, Token *token)
{
  pll_rnode_t *current_node = p->nodes[p->current_node_index];
  assert(!current_node->label);
  current_node->label = token->str;
  token->str = NULL;
}

void destroy_token(Token *token)
{
  free(token->str);
  token->str = NULL;
}

void parse(RTreeParser *p)
{
  Token token;
  Token tokenBL;
  token.str = NULL;
  tokenBL.str = NULL;
  bool end = false;
  rtree_add_node_down(p); // add root
  while (read_token(p, &token)) 
  {
    assert(!end);
    switch(token.type) {
    case TT_SEMICOLON:
      end = true;
      break;
    case TT_COLON:
      read_token(p, &tokenBL);
      assert(tokenBL.type == TT_DOUBLE);
      p->nodes[p->current_node_index]->length = tokenBL.numeric_value;
      destroy_token(&tokenBL);
      break;
    case TT_LEFT_PAR:
      // go down in the tree
      rtree_add_node_down(p);
      break;
    case TT_RIGHT_PAR:
      rtree_go_up(p); 
      break;
    case TT_COMMA:
      rtree_add_node_neighbor(p);
      break;
    case TT_STRING:
    case TT_DOUBLE:
      // todo for now, we allow double label
      rtree_add_label(p, &token); 
      break;
    }
    destroy_token(&token);
  }
  assert(p->parent_node_index == INVALID_INDEX);
}

pll_rtree_t *build_rtree(RTreeParser *p)
{
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
    if (!node->left) {
      node->node_index = tips_index++;
    } else {
      node->node_index = internal_index++;
    }
    tree->nodes[node->node_index] = node;
  }
  tree->root = p->nodes[0];
  tree->root->node_index = internal_index++;
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
  // INIT
  RTreeParser p;
  p.error = error;
  error->type = PET_NOERROR;
  // todo: use const correctly!!
  p.is_file = is_file;
  if (is_file) {
    p.input = getFileContent(input); 
    if (p.input == NULL) {
      error->type = PET_FILE_DO_NOT_EXISTS;
      return NULL;
    }
  } else {
    p.input = (char*)input;
  }
  p.input_current = (char*)p.input;
  p.nodes_capacity = 1000; // todo make dynamic
  p.nodes_number = 0;
  p.nodes = (pll_rnode_t **)
    calloc(p.nodes_capacity, sizeof(pll_rnode_t *));
  assert(p.nodes);
  p.current_node_index = (unsigned int)-1;
  p.depth = 0;
  parse(&p); 
  pll_rtree_t *rtree = build_rtree(&p);
  destroy_rtree_parser(&p);
  if (is_file) {
    free(p.input);
  }
  return rtree;
}


