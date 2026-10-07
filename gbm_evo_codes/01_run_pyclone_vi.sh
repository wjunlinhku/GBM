#!/usr/bin/env bash
set -Eeuo pipefail

# Stage 1: fit PyClone-VI models for all patient-level input files.
# Environment and path configuration.
# PyClone-VI model parameters.
# Input tables are assumed to be prevalidated PyClone-VI input files.

CONDA_SH="${CONDA_SH:-${HOME}/miniconda3/etc/profile.d/conda.sh}"
PYCLONE_ENV="${PYCLONE_ENV:-pyclone-vi}"
[[ -f "$CONDA_SH" ]] || { echo "ERROR: conda setup script not found: $CONDA_SH" >&2; exit 1; }
source "$CONDA_SH"
conda activate "$PYCLONE_ENV"
command -v pyclone-vi >/dev/null 2>&1 || { echo "ERROR: pyclone-vi is unavailable" >&2; exit 1; }

input_dir="${PYCLONE_INPUT_DIR:-/path/to/project/pyclone_input}"
output_dir="${PYCLONE_OUTPUT_DIR:-/path/to/project/pyclone_output}"

mkdir -p "$output_dir"

num_clusters=20
density="beta-binomial"
num_restarts=10

processed_count=0
error_count=0

mapfile -t available_files < <(find "$input_dir" -maxdepth 1 -type f -name "*_pyclone.vi_input.tsv" -print | sort)
if [ ${#available_files[@]} -eq 0 ]; then
    echo "ERROR: No PyClone-VI input files were found in: $input_dir" >&2
    exit 1
fi

echo "Processing ${#available_files[@]} PyClone-VI input files."

for input_file in "${available_files[@]}"; do
    if [ ! -f "$input_file" ]; then
        continue
    fi
    
    patient_id=$(basename "$input_file" "_pyclone.vi_input.tsv")
    h5_output="${output_dir}/${patient_id}_results.h5"
    tsv_output="${output_dir}/${patient_id}_results.tsv"
    
    # Input files are passed directly to PyClone-VI without modification.
    echo "Processing patient: $patient_id"
    if ! pyclone-vi fit -i "$input_file" -o "$h5_output" -c "$num_clusters" -d "$density" -r "$num_restarts"; then
        echo "ERROR: PyClone-VI fit failed for patient: $patient_id" >&2
        error_count=$((error_count + 1)); continue
    fi
    
    if [ ! -f "$h5_output" ]; then
        echo "ERROR: PyClone-VI fit did not produce: $h5_output" >&2
        error_count=$((error_count + 1))
        continue
    fi
    
    if pyclone-vi write-results-file -i "$h5_output" -o "$tsv_output" && [ -f "$tsv_output" ]; then
        processed_count=$((processed_count + 1))
    else
        echo "ERROR: Failed to write results for patient: $patient_id" >&2
        error_count=$((error_count + 1))
    fi
done

echo "Completed PyClone-VI: $processed_count succeeded, $error_count failed. Results: $output_dir"
if [ "$error_count" -gt 0 ]; then exit 1; fi
