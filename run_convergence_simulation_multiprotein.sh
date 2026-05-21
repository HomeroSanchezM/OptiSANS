#!/bin/bash
# Run multiple genetic algorithm simulations for one or more proteins,
# each with the same set of seeds, and generate fitness evolution plots.
# Also generates D2O% vs %D scatter plots for each generation.

set -e  # exit on error

# ---- Timer utilities ----
format_duration() {
    local secs=$1
    printf "%02dh %02dm %02ds" $((secs/3600)) $(( (secs%3600)/60 )) $((secs%60))
}
SCRIPT_START=$(date +%s)

# ---- Fixed parameters (can be changed here) ----
BATCH_SCRIPT="./parallel_process_pdb.sh"
PROCESSES=150
GENERATIONS=500
ELITISM=10
D2O_VAR=5
RATIO_THRESH=0.01

# Seeds to run (42 + 1..9)
SEEDS=(42)

# ---- Reference options ----
# Set NO_DEFAULT_REF=true to skip the automatic protonated-in-D2O / H2O references.
NO_DEFAULT_REF=false

# Add extra reference PDB paths here (space-separated), or leave empty.
# These are passed via --ref to generate_deuterated_pdbs.py.
# Example: REF_PDBS=("original/my_custom_ref1.pdb" "original/my_custom_ref2.pdb")
REF_PDBS=()

# ------------------------------------------------

# Function to display usage
usage() {
    echo "Usage: $0 [--config config.ini] protein1 [protein2 ...]"
    echo ""
    echo "  --config FILE   Path to a config.ini file.  When provided, all GA"
    echo "                  parameters (population, generations, elitism, d2o-var,"
    echo "                  ratio-threshold, restrictions …) are read from the file."
    echo "                  The per-loop --seed is still passed explicitly and"
    echo "                  overrides the INI value, as expected."
    echo "                  When omitted, the hardcoded variables at the top of"
    echo "                  this script are used (legacy behaviour)."
    echo ""
    echo "Examples:"
    echo "  $0 gfp mch rnase"
    echo "  $0 --config config.ini gfp mch"
    echo ""
    echo "Each protein must have a corresponding PDB file at the given path."
    exit 1
}

# ---- Parse optional --config flag and collect protein positional args ----
CONFIG_INI=""
POSITIONAL_ARGS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --config)
            if [[ -z "${2:-}" || "$2" == --* ]]; then
                echo "ERROR: --config requires a file path argument."
                usage
            fi
            CONFIG_INI="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            POSITIONAL_ARGS+=("$1")
            shift
            ;;
    esac
done
# Restore positional parameters to the protein list
set -- "${POSITIONAL_ARGS[@]}"

# Validate --config path if provided
if [ -n "$CONFIG_INI" ]; then
    if [ ! -f "$CONFIG_INI" ]; then
        echo "ERROR: Config file not found: $CONFIG_INI"
        exit 1
    fi
    echo "Using config file: $CONFIG_INI"
fi

# Check if at least one protein is provided
if [ $# -eq 0 ]; then
    echo "ERROR: No protein names given."
    usage
fi

# Build --no_default_ref and --ref arguments to pass to Python
ref_args=()
if [ "$NO_DEFAULT_REF" = true ]; then
    ref_args+=("--no_default_ref")
fi
if [ ${#REF_PDBS[@]} -gt 0 ]; then
    ref_args+=("--ref" "${REF_PDBS[@]}")
fi

for PROTEIN in "$@"; do
    # Extract just the protein name (no directory, no extension)
    PROTEIN_NAME=$(basename "${PROTEIN%.*}")

    echo ""
    echo "========================================="
    echo "Processing protein: $PROTEIN_NAME (from $PROTEIN)"
    echo "========================================="

    INPUT_PDB="${PROTEIN}"
    if [ ! -f "$INPUT_PDB" ]; then
        echo "WARNING: PDB file $INPUT_PDB not found. Skipping $PROTEIN."
        continue
    fi

    BASE_DIR="Rayan_gfp_conc_2.5_gamma_0_tour4_${PROTEIN_NAME}"

    for SEED in "${SEEDS[@]}"; do
        echo ""
        echo "--- Running seed $SEED for $PROTEIN ---"
        SEED_START=$(date +%s)

        # Define output directory for this seed
        OUT_DIR="${BASE_DIR}/seed_${SEED}/${PROTEIN_NAME}/"
        mkdir -p "$OUT_DIR"

        # Build the python command arguments depending on whether a config file
        # was provided.  --batch_script is always passed from the shell script
        # because it is an infrastructure concern, not a GA parameter.
        # --seed is always passed explicitly so that the per-loop value
        # overrides whatever is written in the INI file.
        if [ -n "$CONFIG_INI" ]; then
            # Config-file mode: GA parameters come from the INI; only
            # infrastructure / per-run overrides are added on the CLI.
            python src/python_project/generate_deuterated_pdbs.py \
                "$INPUT_PDB" \
                --config "$CONFIG_INI" \
                --batch_script "$BATCH_SCRIPT" \
                --seed "$SEED" \
                --output_dir "$OUT_DIR" \
                "${ref_args[@]}"
        else
            # Legacy mode: all GA parameters are taken from the variables
            # defined at the top of this script.
            python src/python_project/generate_deuterated_pdbs.py \
                "$INPUT_PDB" \
                --batch_script "$BATCH_SCRIPT" \
                -p "$PROCESSES" \
                -g "$GENERATIONS" \
                -e "$ELITISM" \
                --d2o-var "$D2O_VAR" \
                --seed "$SEED" \
                --ratio-threshold "$RATIO_THRESH" \
                --output_dir "$OUT_DIR" \
                "${ref_args[@]}"
        fi

        # Generate and save the fitness evolution plot
        CSV_FILE="${OUT_DIR}/best_fitness_summary.csv"
        PLOT_FILE="${OUT_DIR}/Fitness_evolution.png"

        if [ -f "$CSV_FILE" ]; then
            python src/python_project/plot_fitness_evolution.py \
                "$CSV_FILE" \
                --annotate \
                --min \
                -o "$PLOT_FILE"
            echo "Plot saved to $PLOT_FILE"
        else
            echo "Warning: $CSV_FILE not found – skipping plot."
        fi

        # Generate D2O% vs %D scatter plots for each generation
        python src/python_project/d2o_vs_d.py "$OUT_DIR"
        echo "D2O vs %D plots saved in ${OUT_DIR}/generation_plots_d2o_vs_d/"

        SEED_END=$(date +%s)
        echo "Seed $SEED elapsed: $(format_duration $((SEED_END - SEED_START)))"
    done
done

echo ""
echo "All simulations completed."
SCRIPT_END=$(date +%s)
echo "Total elapsed time: $(format_duration $((SCRIPT_END - SCRIPT_START)))"
