#!/usr/bin/env bash
set -Eeuo pipefail

# Stage 6: run RECAP across the requested range of model complexities.

# RECAP path and model configuration.
INPUT_FILE="${RECAP_INPUT_FILE:-/path/to/project/recap/recap_input.txt}"
OUTPUT_BASE="${RECAP_OUTPUT_DIR:-/path/to/project/recap/recap_output}"
RECAP_BINARY="${RECAP_BINARY:-/path/to/RECAP/build/recap}"
ROOT_GENE="GL"
RESTART_COUNT=5
K_MIN=2
K_MAX=6
CUTOFF=0
if [ "${CONDA_DEFAULT_ENV:-}" != "recap" ]; then
    echo "ERROR: Activate the recap conda environment before running this script." >&2
    exit 1
fi

if [ ! -f "$INPUT_FILE" ] || [ ! -f "$RECAP_BINARY" ]; then
    echo "ERROR: RECAP input file or binary was not found." >&2
    exit 1
fi

if ! command -v dot &> /dev/null; then
    echo "ERROR: Graphviz 'dot' is not available on PATH." >&2
    exit 1
fi

mkdir -p "$OUTPUT_BASE"

successful_runs=()
failed_runs=()

run_recap_k() {
    local k=$1
    local output_prefix="${OUTPUT_BASE}/cosmic_k${k}"
    local log_file="${OUTPUT_BASE}/recap_k${k}_log.txt"
    
    timeout $TIMEOUT_SECONDS $RECAP_BINARY \
        -k $k -p "$output_prefix" -R "$ROOT_GENE" \
        "$INPUT_FILE" -r $RESTART_COUNT -c $CUTOFF > "$log_file" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "RECAP completed for k=${k}."
        successful_runs+=($k)
    else
        echo "ERROR: RECAP failed for k=${k}; see $log_file" >&2
        failed_runs+=($k)
    fi
}

for k in $(seq $K_MIN $K_MAX); do
    run_recap_k $k
done

if [ ${#successful_runs[@]} -gt 0 ]; then
    cd "$OUTPUT_BASE" || exit
    dot_files=$(ls *.dot 2>/dev/null)
    
    if [ -z "$dot_files" ]; then
        echo "WARNING: No Graphviz .dot files were produced." >&2
    else
        for dot_file in $dot_files; do
            base_name="${dot_file%.dot}"
            dot -Tpng "$dot_file" -o "${base_name}.png"
            dot -Tpdf "$dot_file" -o "${base_name}.pdf"
            dot -Tsvg "$dot_file" -o "${base_name}.svg"
        done
    fi

    cat > analyze_results.py << 'EOF'
import os, sys, pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from datetime import datetime

def main():
    output_base = sys.argv[1]
    k_list = [int(x) for x in sys.argv[2].split(',') if x]
    results = []

    for k in k_list:
        sum_f = os.path.join(output_base, f"cosmic_k{k}_summary.tsv")
        if os.path.exists(sum_f):
            df = pd.read_csv(sum_f, sep='\t')
            if not df.empty:
                row = df.iloc[0]
                res = {'k': k, 'cost': row.get('Cost', 0), 'time': row.get('Time', 0)}
                for i in range(k):
                    res[f'c{i}'] = row.get(f'Cluster {i}', 0)
                results.append(res)
    
    if not results:
        print("Error: No valid data found.")
        return
    
    rdf = pd.DataFrame(results).sort_values('k')
    
    rdf['diff'] = rdf['cost'].shift(1) - rdf['cost']
    rdf['rate'] = rdf['diff'] / rdf['cost'].shift(1)
    
    optimal_k = rdf['k'].iloc[0]
    for i in range(1, len(rdf)):
        if rdf.iloc[i]['rate'] < 0.05:
            optimal_k = rdf.iloc[i-1]['k']
            break
    if not optimal_k: optimal_k = rdf['k'].max()

    plt.figure(figsize=(10, 6))
    plt.plot(rdf['k'], rdf['cost'], 'o-', color='blue', lw=2)
    plt.axvline(x=optimal_k, color='red', linestyle='--', label=f'Best k={optimal_k}')
    plt.title('Cost vs K-clusters')
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.savefig('recap_analysis_plot.png', dpi=300)
    
    with open('RECAP_Report.md', 'w') as f:
        f.write(f"# RECAP \n\n: {datetime.now()}\n\n k : **{optimal_k}**\n\n")
        f.write(rdf.to_markdown(index=False))

    print(optimal_k)

if __name__ == "__main__":
    main()
EOF

    optimal_k=$(python3 analyze_results.py "$OUTPUT_BASE" "$(IFS=,; echo "${successful_runs[*]}")")

    echo "RECAP completed. Optimal k: $optimal_k. Results: $OUTPUT_BASE"
else
    echo "ERROR: RECAP failed for all requested k values." >&2
    exit 1
fi
