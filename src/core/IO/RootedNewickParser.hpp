#pragma once

#include <trees/PLLRootedTree.hpp>

#include <IO/NewickParserCommon.hpp>
/**
 *  Parse a rooted tree from a newick string or file
 *  If an error occures, returns NULL and fills the 
 *  error object.
 */
pll_rtree_t * custom_rtree_parse_newick(const char *s, 
    bool is_file,
    ParsingError *error);




