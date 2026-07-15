#!/bin/bash

# Usage:
#   Full mode (generation 0): process ALL PDB files in the folder + ref
#       $0 <deuterated_pdbs_folder> [num_jobs] [--d2o X.XX]
#
#   Incremental mode (generation N > 0): process only the listed files + ref
#       $0 <deuterated_pdbs_folder> <pdb_list_file> [num_jobs] [--d2o X.XX]
#
#   <pdb_list_file> is a plain text file with one absolute PDB path per line.
#   When it is supplied the ref/ subfolder is NOT re-processed (outputs already exist).
#
# D2O flag selection:
#   --d2o X.XX   use this explicit value for every PDB file
#   otherwise    extract D2O from the filename pattern:
#                *_d2o<digits> -> --d2o <digits>/100
#
# Reference PDB naming conventions for inferred D2O mode:
#   *_total_deuteration  -> --d2o 1   (protonated in D2O)
#   *_total_protonation  -> --d2o 0   (protonated in H2O)
#   *_d2oXX              -> --d2o XX/100  (any other ref with explicit D2O value)
#   otherwise            -> no --d2o flag

# ---------------------------------------------------------------------------
# Check mandatory argument
# ---------------------------------------------------------------------------
if [ $# -eq 0 ]; then
    echo "Usage: $0 <deuterated_pdbs_folder> [pdb_list_file] [num_jobs] [--d2o X.XX]"
    exit 1
fi

SECONDS=0

input_dir="$1"
pdb_list_file=""
num_jobs=150   # default parallel jobs
explicit_d2o=""
conc=2.5       # default concentration
export explicit_d2o conc

# ---------------------------------------------------------------------------
# Parse optional arguments
# ---------------------------------------------------------------------------
while [ $# -gt 1 ]; do
    case "$2" in
        --d2o)
            if [ -n "$3" ] && [ "$3" != --* ]; then
                explicit_d2o="$3"
                export explicit_d2o
                shift 2
                continue
            else
                echo "Error: --d2o requires a numeric value" >&2
                exit 1
            fi
            ;;
        --conc)
            if [ -n "$3" ] && [ "$3" != --* ]; then
                conc="$3"
                export conc
                shift 2
                continue
            else
                echo "Error: --conc requires a numeric value" >&2
                exit 1
            fi
            ;;
        *)
            if [ -z "$pdb_list_file" ] && [ -f "$2" ]; then
                pdb_list_file="$2"
            elif [ -z "$pdb_list_file" ] && [[ "$2" =~ ^[0-9]+$ ]]; then
                num_jobs="$2"
            fi
            ;;
    esac
    shift
done

# ---------------------------------------------------------------------------
# Validate input directory
# ---------------------------------------------------------------------------
if [ ! -d "$input_dir" ]; then
    echo "Error: The folder '$input_dir' does not exist!"
    exit 1
fi

# ---------------------------------------------------------------------------
# Compute output directory for SANS results (primus_out)
#   - place it in the same parent as the input folder
#   - if the folder name ends with '_deuterated_pdbs', strip that suffix
#     before appending '_primus_out'; otherwise use the full folder name.
# ---------------------------------------------------------------------------
parent_dir=$(dirname "$input_dir")
folder_name=$(basename "$input_dir")
prefix="${folder_name%_deuterated_pdbs}"
output_dir="${parent_dir}/${prefix}_primus_out"

# Create or clean output directories
if [ -d "$output_dir" ] && [ -z "$pdb_list_file" ]; then
    # Full mode: wipe the output directory (except ref/)
    echo "=== Full mode: cleaning output directory ==="
    find "$output_dir" -maxdepth 1 -name "*.dat" -delete
    find "$output_dir" -maxdepth 1 -name "*.out" -delete
    echo ""
fi

mkdir -p "$output_dir"
mkdir -p "$output_dir/ref"

echo "=== Configuration ==="
echo "Input folder:    $input_dir"
echo "Output folder:   $output_dir"
echo "Prefix:          $prefix"
echo "Parallel jobs:   $num_jobs"
echo "Concentration:   $conc"
if [ -n "$explicit_d2o" ]; then
    echo "D2O override:    $explicit_d2o"
fi
if [ -n "$pdb_list_file" ]; then
    echo "Mode:            INCREMENTAL (file list: $pdb_list_file)"
else
    echo "Mode:            FULL (all PDB files)"
fi
echo ""

# ---------------------------------------------------------------------------
# Determine the --d2o flag for a reference PDB file.
# Priority order:
#   1) explicit CLI value --d2o X.XX
#   2) naming convention:
#        *_total_deuteration  -> --d2o 1
#        *_total_protonation  -> --d2o 0
#        *_d2oXX              -> --d2o XX/100
#   3) otherwise no flag
# ---------------------------------------------------------------------------
_get_ref_d2o_flag() {
    local basename="$1"
    if [ -n "$explicit_d2o" ]; then
        echo "--d2o $explicit_d2o"
        return
    fi
    if [[ "$basename" == *"_total_deuteration" ]]; then
        echo "--d2o 1"
    elif [[ "$basename" == *"_total_protonation" ]]; then
        echo "--d2o 0"
    elif [[ "$basename" =~ _d2o([0-9]+) ]]; then
        local d2o_int=$((10#${BASH_REMATCH[1]}))
        local d2o_value
        d2o_value=$(LC_NUMERIC=C awk "BEGIN { printf \"%.2f\", $d2o_int / 100 }")
        echo "--d2o $d2o_value"
    else
        echo ""
    fi
}
export -f _get_ref_d2o_flag

process_ref_file() {
    pdb_file="$1"
    output_dir="$2"

    basename=$(basename "$pdb_file" .pdb)
    output_file="$output_dir/ref/${basename}.dat"

    d2o_flag=$(_get_ref_d2o_flag "$basename")

    echo "Processing ref: $pdb_file -> $output_file (d2o flag: '${d2o_flag:-none}')"
    # shellcheck disable=SC2086
    ./Pepsi-SANS-Linux/Pepsi-SANS "$pdb_file" --hModel 3 --conc "$conc" $d2o_flag -o "$output_file"
}
export -f process_ref_file

# ---------------------------------------------------------------------------
# process one MAIN file
# ---------------------------------------------------------------------------
process_main_file() {
    pdb_file="$1"
    output_dir="$2"

    basename=$(basename "$pdb_file" .pdb)
    output_file="$output_dir/${basename}.dat"

    # Use explicit --d2o if provided; otherwise infer from filename.
    if [ -n "$explicit_d2o" ]; then
        d2o_flag="--d2o $explicit_d2o"
    elif [[ "$basename" =~ _d2o([0-9]+) ]]; then
        d2o_int=$((10#${BASH_REMATCH[1]}))
        d2o_value=$(LC_NUMERIC=C awk "BEGIN { printf \"%.2f\", $d2o_int / 100 }")
        d2o_flag="--d2o $d2o_value"
    else
        echo "Warning: could not extract d2o value from '$basename', skipping --d2o flag"
        d2o_flag=""
    fi

    echo "Processing: $pdb_file -> $output_file (d2o flag: '${d2o_flag:-none}')"
    # shellcheck disable=SC2086
    ./Pepsi-SANS-Linux/Pepsi-SANS "$pdb_file" --hModel 3 --conc "$conc" $d2o_flag -o "$output_file"
}
export -f process_main_file

# ---------------------------------------------------------------------------
# Process ref/ subfolder — only in full mode
# ---------------------------------------------------------------------------
if [ -z "$pdb_list_file" ]; then
    if [ -d "$input_dir/ref" ]; then
        ref_files=("$input_dir/ref"/*.pdb)
        if [ -e "${ref_files[0]}" ]; then
            echo "=== Processing ${#ref_files[@]} reference file(s) (parallel -j $num_jobs) ==="
            printf '%s\n' "${ref_files[@]}" | parallel -j "$num_jobs" -k \
                process_ref_file {} "$output_dir"
            echo ""
        else
            echo "No .pdb files found in $input_dir/ref/"
            echo ""
        fi
    fi
else
    echo "=== Incremental mode: skipping ref/ (already processed) ==="
    echo ""
fi

# ---------------------------------------------------------------------------
# Process main files
# ---------------------------------------------------------------------------
if [ -n "$pdb_list_file" ]; then
    # Incremental mode: read the list of new PDB files
    mapfile -t main_files < "$pdb_list_file"
    # Filter to only existing files
    existing_files=()
    for f in "${main_files[@]}"; do
        [ -f "$f" ] && existing_files+=("$f")
    done
    if [ ${#existing_files[@]} -eq 0 ]; then
        echo "No valid PDB files found in list file '$pdb_list_file'"
    else
        echo "=== Processing ${#existing_files[@]} new file(s) (parallel -j $num_jobs) ==="
        printf '%s\n' "${existing_files[@]}" | parallel -j "$num_jobs" -k \
            process_main_file {} "$output_dir"
    fi
else
    # Full mode: process all PDB files in input_dir
    main_files=("$input_dir"/*.pdb)
    if [ -e "${main_files[0]}" ]; then
        echo "=== Processing ${#main_files[@]} main file(s) (parallel -j $num_jobs) ==="
        printf '%s\n' "${main_files[@]}" | parallel -j "$num_jobs" -k \
            process_main_file {} "$output_dir"
    else
        echo "No .pdb files found in $input_dir/"
    fi
fi

echo ""
echo "=== Simulation complete! ==="
echo "Results saved in: $output_dir"

elapsed=$SECONDS
hours=$((elapsed / 3600))
minutes=$(( (elapsed % 3600) / 60 ))
seconds=$((elapsed % 60))
printf "Total execution time: %02dh %02dm %02ds\n" "$hours" "$minutes" "$seconds"
