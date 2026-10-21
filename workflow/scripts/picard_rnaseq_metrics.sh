#!/bin/bash
#
# Wrapper for Picard CollectRnaSeqMetrics with optional multi-threading.
#
# Usage: picard_rnaseq_metrics.sh <threads> <CollectRnaSeqMetrics arguments>
#
# threads <= 1: runs 'picard CollectRnaSeqMetrics' unchanged.
# threads  > 1: splits the (coordinate-sorted, indexed) input BAM into genomic chunks,
#               processes them in parallel and merges the results into the same
#               metrics file (see ParallelRnaSeqMetrics.java).
#
# Must run inside the Picard container. Java heap is set with _JAVA_OPTIONS, e.g. -Xmx25g.

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo ""
  echo "Picard CollectRnaSeqMetrics, optionally multi-threaded."
  echo ""
  echo "Usage: $0 <threads> <CollectRnaSeqMetrics arguments>"
  echo ""
  exit 1
fi

THREADS=$1
shift

PICARD=$(command -v picard)
if [[ $THREADS -le 1 ]]; then
    exec "$PICARD" CollectRnaSeqMetrics "$@"
fi

# Same JAR, Java and locale as the 'picard' launcher script
export LC_ALL=en_US.UTF-8
JAR_DIR=$(dirname "$(readlink -f "$PICARD")")
JAVA=$(dirname "$(dirname "$JAR_DIR")")/bin/java
[[ -x $JAVA ]] || JAVA=java
SCRIPT_DIR=$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")

exec "$JAVA" -Dthreads="$THREADS" -cp "$JAR_DIR/picard.jar" "$SCRIPT_DIR/ParallelRnaSeqMetrics.java" "$@"
