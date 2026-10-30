#!/usr/bin/env bash
# Pack the genomes for Ridom SeqSphere+ (scripts/01): the pool (175 complex isolates incl. M055) + 7 GTDB references, as tar.gz and zip.
#   sbatch submit_export.sbatch
# Inputs : data/imm/screening_map_results.csv, runs/mkan329/assembly/results, GTDB_REPRESENTATIVES
# Outputs: output/kansasii-ridom/ridom_assemblies.{tar.gz,zip} (+ ridom_assemblies.MANIFEST.txt)
#SBATCH --job-name=kansasii-ridom-export
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/shares/sander.imm.uzh/MM/kansasii/output/kansasii-ridom/logs/%x-%j.log

set -euo pipefail
source "${KANSASII_ROOT:-/shares/sander.imm.uzh/MM/kansasii}/repos/kansasii-docs/templates/kansasii_paths.sh"
REPO=$REPOS/kansasii-ridom
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate env_ridom
cd "$REPO/scripts"
python 01_export_ridom.py
echo "done $(date +%T)"

# provenance (see kansasii-docs/docs/PROVENANCE.md); loose files in the output root get <file>.MANIFEST.txt
MANIFEST_NAME=ridom_assemblies.MANIFEST.txt bash "$MANIFEST_SCRIPT" "$OUT_RIDOM" "$REPO" submit_export.sbatch \
  "$SCREENING_MAP_RESULTS" "$ASSEMBLIES" "$GTDB_REPRESENTATIVES"
