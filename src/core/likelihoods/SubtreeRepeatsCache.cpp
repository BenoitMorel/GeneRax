#include "SubtreeRepeatsCache.hpp"
#include <likelihoods/reconciliation_models/AbstractReconciliationModel.hpp>

static unsigned int integerHash(unsigned int x) {
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = ((x >> 16) ^ x) * 0x45d9f3b;
  x = (x >> 16) ^ x;
  return x;
}

unsigned long hashing_func_subtree::operator()(const pll_unode_t *subtree) const {
  if (!subtree->next) {
    return static_cast<unsigned long>(_cache.geneToSpecies(subtree)); 
  }
  auto left = subtree->next->back;
  auto right = subtree->next->next->back;
  auto leftHash = integerHash(_cache.getRepeatIndexNoCheck(left));
  auto rightHash = integerHash(_cache.getRepeatIndexNoCheck(right));
  return static_cast<unsigned long>(integerHash(leftHash + rightHash));
}

bool key_equal_fn_subtree::operator()(const pll_unode_t *subtree1, const pll_unode_t *subtree2) const {
  if (!subtree1->next && !subtree2->next) {
    return _cache.geneToSpecies(subtree1) == _cache.geneToSpecies(subtree2);
  } else if (!subtree1->next || !subtree2->next) {
    return false;
  }
  auto rid1_left = _cache.getRepeatIndexNoCheck(subtree1->next->back);
  auto rid1_right = _cache.getRepeatIndexNoCheck(subtree1->next->next->back);
  auto rid2_left = _cache.getRepeatIndexNoCheck(subtree2->next->back);
  auto rid2_right = _cache.getRepeatIndexNoCheck(subtree2->next->next->back);
  bool res1 = (rid1_left == rid2_left) && (rid1_right == rid2_right);
  bool res2 = (rid1_left == rid2_right) && (rid1_right == rid2_left);
  return res1 || res2;
}

