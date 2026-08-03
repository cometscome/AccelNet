#!/usr/bin/env bash
set -euo pipefail

linux_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
result_root=${BENCH_RESULT_DIR:-$linux_dir/results/variant-run-$(date +%Y%m%d-%H%M%S)}
mkdir -p "$result_root"
printf 'STARTING\n' > "$result_root/status.txt"

if command -v setsid >/dev/null 2>&1; then
    nohup setsid env BENCH_RESULT_DIR="$result_root" \
        "$linux_dir/run_variant_benchmarks.sh" \
        > "$result_root/driver.log" 2>&1 < /dev/null &
else
    nohup env BENCH_RESULT_DIR="$result_root" \
        "$linux_dir/run_variant_benchmarks.sh" \
        > "$result_root/driver.log" 2>&1 < /dev/null &
fi
background_pid=$!
printf '%s\n' "$background_pid" > "$result_root/background.pid"

echo "Background benchmark started."
echo "Result directory: $result_root"
echo "PID: $background_pid"
echo "Status: cat '$result_root/status.txt'"
echo "Progress: tail -f '$result_root/driver.log'"
