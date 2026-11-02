#!/usr/bin/env bash
# Check the Ridom exports in data/imm (scripts/00) and derive all tables for the paper from them (02-08; the archive for Ridom is submit_export.sbatch).
#   sbatch submit_all.sbatch
# Inputs : data/imm/ridom-*, data/imm/screening_map_results.csv
# Outputs: output/kansasii-ridom/results/ export_status, distance_*, scheme_*, qc*, locus_stats, target_failures, cluster*, mst_*, tree_upgma_*, control_distances, paper_*
#SBATCH --job-name=kansasii-ridom
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=00:20:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/kansasii-ridom/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/kansasii-ridom
OUTPUT=$OUT_RIDOM/results
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_ridom
cd "$REPO/scripts"
# A failing check (exit 1) must not hide the other tables: run every step, remember failures, exit 1 at the end.
# 08 reads the results of 02-07, so it only runs when all of them succeeded (no stale files).
rc=0
run() { echo "=== $1.py ($(date +%T))"; python "$1.py" || { echo "FAILED: $1.py"; rc=1; }; }
run 00_check_exports
for s in 02_distance_matrices 03_scheme 04_qc 05_clusters 06_mst_trees 07_control_distances; do run $s; done
if [ "$rc" -eq 0 ]; then run 08_paper_numbers; else echo "SKIPPED: 08_paper_numbers.py (earlier step failed)"; fi
echo "done $(date +%T) (exit $rc)"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT" "$REPO" submit_all.sbatch "$SCREENING_MAP_RESULTS" "$IMM_DIR"
exit $rc
