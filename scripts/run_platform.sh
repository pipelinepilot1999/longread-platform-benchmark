#!/usr/bin/env bash
# Launch (or resume) one platform run detached from the SSH session, pinning the Nextflow session id.
#   scripts/run_platform.sh ont  data/ont_chr20.raw.bam  [--target_cov 49.61] [extra nextflow args]
#   scripts/run_platform.sh hifi data/hifi_chr20.raw.bam --target_cov 49.61
#   RESUME=<session-uuid> scripts/run_platform.sh ont ... -> resumes that exact session
set -euo pipefail
platform=$1; raw_bam=$2; shift 2
cd "$(dirname "$0")/.."
mkdir -p logs
stamp=$(date +%Y%m%d_%H%M%S)
log=logs/nf_${platform}_${stamp}.log
resume_arg=""
if [ -n "${RESUME:-}" ]; then resume_arg="-resume ${RESUME}"; fi
cmd="nextflow -log logs/.nextflow_${platform}_${stamp}.log run main.nf --platform ${platform} --raw_bam ${raw_bam} $* ${resume_arg} -with-trace -ansi-log false"
echo "$cmd" | tee -a "$log"
nohup setsid bash -c "$cmd" >> "$log" 2>&1 < /dev/null &
sleep 12
sid=$(grep -m1 -o "Session UUID: [0-9a-f-]*" "logs/.nextflow_${platform}_${stamp}.log" 2>/dev/null | awk '{print $3}')
echo "session id (use RESUME=<id> to resume this run): ${sid:-not yet in log; grep 'Session UUID' logs/.nextflow_${platform}_${stamp}.log}" | tee -a "$log"
echo "log: $log"
