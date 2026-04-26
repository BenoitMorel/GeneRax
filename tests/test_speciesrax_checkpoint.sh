#!/bin/bash
# Test script for SpeciesRax checkpointing
# Uses 100 randomly selected families from cyano_simulated dataset
#
# Verifies that:
# 1. A fresh run (no --checkpoint) produces identical results to an
#    uninterrupted --checkpoint run
# 2. An interrupted (phase-level) + resumed --checkpoint run produces
#    identical results to a fresh run
# 3. An interrupted (iteration-level) + resumed --checkpoint run
#    completes successfully (result may differ slightly because the
#    evaluator's internal state cannot be fully serialized)
#
# Usage: Run from the GeneRax repository root:
#   bash tests/test_speciesrax_checkpoint.sh

set -euo pipefail

GENERAX=build/bin/generax
DATADIR=generax_data/cyano_simulated
TESTDIR=tests/speciesrax_checkpoint_test
FAMILIES_FILE=$TESTDIR/families_100.txt

if [ ! -f "$GENERAX" ]; then
  echo "Error: Cannot find $GENERAX. Please build GeneRax first (./install.sh)."
  exit 1
fi

if [ ! -d "$DATADIR" ]; then
  echo "Error: Cannot find $DATADIR. Please download and extract the test data first:"
  echo "  wget https://cme.h-its.org/exelixis/material/generax_data.tar.gz"
  echo "  tar xzf generax_data.tar.gz"
  exit 1
fi

# Create test directory
mkdir -p "$TESTDIR"

# Generate families file with 100 randomly selected families
echo "Generating families file with 100 random families..."
python3 -c "
import os, random
random.seed(42)
families_dir = '$DATADIR/families'
families = sorted(os.listdir(families_dir))
selected = sorted(random.sample(families, 100))
print('[FAMILIES]')
for fam in selected:
    gene_tree = os.path.join(families_dir, fam, 'gene_trees', 'raxml-ng.LG+G+I.geneTree.newick')
    mapping = os.path.join(families_dir, fam, 'mappings', 'mapping.link')
    if os.path.exists(gene_tree) and os.path.exists(mapping):
        print(f'- {fam}')
        print(f'starting_gene_tree = {gene_tree}')
        print(f'mapping = {mapping}')
" > "$FAMILIES_FILE"

NUM_FAMILIES=$(grep -c "^- " "$FAMILIES_FILE")
echo "Generated families file with $NUM_FAMILIES families"

COMMON_ARGS="--families $FAMILIES_FILE --species-tree MiniNJ --strategy SKIP --rec-model UndatedDTL --per-family-rates --prune-species-tree --si-strategy HYBRID"

# Clean up any previous test output
rm -rf "$TESTDIR"/output_*

FAILED=0

echo ""
echo "===== Test 1: Fresh run (no --checkpoint) ====="
echo ""

$GENERAX $COMMON_ARGS --prefix "$TESTDIR/output_fresh"

echo ""
echo "===== Test 2: Uninterrupted --checkpoint run ====="
echo ""

$GENERAX $COMMON_ARGS --prefix "$TESTDIR/output_checkpoint" --checkpoint

# Compare fresh vs uninterrupted checkpoint (must be identical)
echo ""
echo "--- Comparing fresh vs uninterrupted checkpoint ---"
if diff -q "$TESTDIR/output_fresh/species_trees/inferred_species_tree.newick" \
           "$TESTDIR/output_checkpoint/species_trees/inferred_species_tree.newick" > /dev/null 2>&1; then
  echo "PASS: Species trees identical"
else
  echo "FAIL: Species trees differ"
  FAILED=1
fi

if diff -q "$TESTDIR/output_fresh/stats.txt" \
           "$TESTDIR/output_checkpoint/stats.txt" > /dev/null 2>&1; then
  echo "PASS: Stats identical"
else
  echo "FAIL: Stats differ"
  FAILED=1
fi

echo ""
echo "===== Test 3: Phase-level interrupt + resume ====="
echo "(kill before iteration checkpoint is saved)"
echo ""

# 10 seconds is enough to save the starting species tree (phase 1)
# but not enough to complete any search iteration
timeout 10 $GENERAX $COMMON_ARGS --prefix "$TESTDIR/output_phase_resume" --checkpoint || true

echo ""
echo "Checkpoint state after interruption:"
cat "$TESTDIR/output_phase_resume/checkpoint.txt" 2>/dev/null || echo "(no checkpoint file)"
echo ""

# Resume
$GENERAX $COMMON_ARGS --prefix "$TESTDIR/output_phase_resume" --checkpoint

echo ""
echo "--- Comparing fresh vs phase-level resume ---"
if diff -q "$TESTDIR/output_fresh/species_trees/inferred_species_tree.newick" \
           "$TESTDIR/output_phase_resume/species_trees/inferred_species_tree.newick" > /dev/null 2>&1; then
  echo "PASS: Species trees identical"
else
  echo "FAIL: Species trees differ"
  FAILED=1
fi

if diff -q "$TESTDIR/output_fresh/stats.txt" \
           "$TESTDIR/output_phase_resume/stats.txt" > /dev/null 2>&1; then
  echo "PASS: Stats identical"
else
  echo "FAIL: Stats differ"
  FAILED=1
fi

echo ""
echo "===== Test 4: Iteration-level interrupt + resume ====="
echo "(kill after at least one iteration checkpoint is saved)"
echo ""

# 25 seconds should be enough for 1-2 iteration checkpoints
timeout 25 $GENERAX $COMMON_ARGS --prefix "$TESTDIR/output_iter_resume" --checkpoint || true

echo ""
echo "Checkpoint state after interruption:"
cat "$TESTDIR/output_iter_resume/checkpoint.txt" 2>/dev/null || echo "(no checkpoint file)"
echo ""

if grep -q "hybrid_iteration" "$TESTDIR/output_iter_resume/checkpoint.txt" 2>/dev/null; then
  echo "Iteration-level checkpoint found, resuming..."
  $GENERAX $COMMON_ARGS --prefix "$TESTDIR/output_iter_resume" --checkpoint

  echo ""
  echo "--- Checking iteration-level resume ---"
  if [ -f "$TESTDIR/output_iter_resume/species_trees/inferred_species_tree.newick" ]; then
    echo "PASS: Completed successfully"
  else
    echo "FAIL: No output species tree"
    FAILED=1
  fi

  if [ -f "$TESTDIR/output_iter_resume/stats.txt" ]; then
    echo "PASS: Stats file exists"
  else
    echo "FAIL: No stats file"
    FAILED=1
  fi
else
  echo "NOTE: Iteration-level checkpoint was not reached (timing-dependent)."
  echo "      Skipping iteration-level resume test."
fi

echo ""
echo "===== Results ====="
echo ""
echo "Reconciliation likelihoods:"
echo "  Fresh:                $(grep RecLL "$TESTDIR/output_fresh/stats.txt")"
echo "  Checkpoint:           $(grep RecLL "$TESTDIR/output_checkpoint/stats.txt")"
echo "  Phase-level resume:   $(grep RecLL "$TESTDIR/output_phase_resume/stats.txt")"
if [ -f "$TESTDIR/output_iter_resume/stats.txt" ]; then
  echo "  Iter-level resume:    $(grep RecLL "$TESTDIR/output_iter_resume/stats.txt")"
fi

# Clean up
rm -rf "$TESTDIR"

echo ""
if [ $FAILED -eq 0 ]; then
  echo "All tests PASSED!"
  exit 0
else
  echo "Some tests FAILED!"
  exit 1
fi
