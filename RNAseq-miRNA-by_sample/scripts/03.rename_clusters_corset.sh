#!/bin/bash

# Renames Corset clusters to miRNA names (POSIX-compatible)
# Usage: ./03.rename_clusters_corset.sh [clusters_file] [counts_file] [output_file]
#
# Example (individual):
#   ./03.rename_clusters_corset.sh results/corset_mature/mature-clusters.txt \
#                                    results/corset_mature/mature-counts.txt \
#                                    results/count_featurecounts/mature-counts-gga.txt

CLUSTERS_FILE=${1:-"mature-clusters.txt"}
COUNTS_FILE=${2:-"mature-counts.txt"}
OUTPUT_FILE=${3:-"mature-counts-gga.txt"}

if [ ! -f "$CLUSTERS_FILE" ] || [ ! -f "$COUNTS_FILE" ]; then
    echo "Error: Input files not found."
    exit 1
fi

awk -F'\t' '
NR==FNR {
    if ($0 == "" || $0 ~ /^\[/) next
    mir_id = $1
    cluster_id = $2
    map[cluster_id] = map[cluster_id] SUBSEP mir_id
    next
}

FNR==1 {
    print $0
    next
}

{
    cluster = $1
    if (cluster in map) {
        n = split(map[cluster], arr, SUBSEP)
        for (i=1; i<=n; i++) {
            printf "%s", arr[i]
            for (j=2; j<=NF; j++) {
                printf "\t%s", $j
            }
            printf "\n"
        }
    } else {
        print $0
    }
}
' "$CLUSTERS_FILE" "$COUNTS_FILE" > "$OUTPUT_FILE"
