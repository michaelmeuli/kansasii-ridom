#!/usr/bin/env bash
# Check the Ridom exports in data/imm against the list of exports still needed (scripts/00 only; the archive for Ridom is submit_export.sbatch).
#   sbatch submit_all.sbatch
# Inputs : data/imm/ridom-*, data/imm/screening_map_results.csv
# Outputs: output/kansasii-ridom/results/export_status.tsv
#SBATCH --job-name=kansasii-ridom
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:10:00
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
echo "done $(date +%T)"

# provenance (see kansasii-docs/docs/PROVENANCE.md)
bash "$MANIFEST_SCRIPT" "$OUTPUT" "$REPO" submit_all.sbatch "$SCREENING_MAP_RESULTS" "$IMM_DIR"
