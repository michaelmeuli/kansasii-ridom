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
echo "=== 00_check_exports.py ($(date +%T))"
python 00_check_exports.py
echo "=== 02_distance_matrices.py ($(date +%T))"
python 02_distance_matrices.py
for s in 03_scheme 04_qc 05_clusters 06_mst_trees 07_control_distances 08_paper_numbers; do
  echo "=== $s.py ($(date +%T))"
  python $s.py
done
echo "done $(date +%T)"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT" "$REPO" submit_all.sbatch "$SCREENING_MAP_RESULTS" "$IMM_DIR"
