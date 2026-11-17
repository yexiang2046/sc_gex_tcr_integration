#!/bin/bash
# =============================================================================
# sync_cellranger.sh
# Rsync 10x Cell Ranger "outs/" directory while skipping large files
# not needed for typical downstream analysis (Seurat, Scanpy, etc.)
#
# Large files skipped:
#   - possorted_genome_bam.bam + .bai (tens of GB)
#   - raw_feature_bc_matrix/ + raw_feature_bc_matrix.h5
#   - molecule_info.h5 (and variants like raw_molecule_info.h5)
#   - Any other *.bam / *.bai files
#
# Files kept (essential for downstream):
#   - filtered_feature_bc_matrix/ + .h5
#   - analysis/ (clustering, UMAP, etc.)
#   - web_summary.html + metrics_summary.csv
#   - cloupe.cloupe (optional but often useful for Loupe Browser)
#
# Usage:
#   ./sync_cellranger.sh [-n] <source> <destination>
#     -n   Dry-run mode (shows what would be transferred, no changes)
#
# Examples:
#   ./sync_cellranger.sh /data/run123/outs/ /local/analysis/run123/
#   ./sync_cellranger.sh -n user@cluster:/data/run123/outs/ /local/analysis/run123/
# =============================================================================

set -euo pipefail

# Parse optional -n flag for dry-run
DRY_RUN=""
while getopts ":n" opt; do
  case "${opt}" in
    n) DRY_RUN="--dry-run" ;;
    \?) echo "ERROR: Invalid option" >&2; exit 1 ;;
  esac
done
shift $((OPTIND - 1))

if [ $# -ne 2 ]; then
  echo "Usage: $0 [-n] <source_dir> <destination_dir>" >&2
  echo "  -n : dry-run (show changes without transferring)" >&2
  exit 1
fi

SOURCE="${1%/}"   # remove trailing slash if present
DEST="${2%/}"

# Basic validation
if [ ! -d "${SOURCE}" ] && [[ "${SOURCE}" != *:* ]]; then
  echo "ERROR: Source directory does not exist: ${SOURCE}" >&2
  exit 1
fi

echo "=== Cell Ranger rsync started ==="
echo "Source      : ${SOURCE}"
echo "Destination : ${DEST}"
[ -n "${DRY_RUN}" ] && echo "MODE        : DRY-RUN (nothing will be transferred)"
echo "Excluding large/unnecessary files for downstream analysis..."
echo

rsync -avh --progress ${DRY_RUN} \
  --exclude='**/possorted_genome_bam.bam*' \
  --exclude='**/raw_feature_bc_matrix/' \
  --exclude='**/raw_feature_bc_matrix.h5' \
  --exclude='**/molecule_info.h5' \
  --exclude='**/raw_molecule_info.h5' \
  --exclude='**/*.bam' \
  --exclude='**/*.bai' \
  "${SOURCE}/" "${DEST}/"

echo
echo "=== Sync completed! ==="
echo "Only files needed for downstream single-cell analysis were transferred."
