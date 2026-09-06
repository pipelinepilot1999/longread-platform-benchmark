#!/usr/bin/env bash
# Disk guard: exit 75 (-> Nextflow 'terminate') if free space on the work volume is below the threshold in GB.
min_gb=${1:-40}
free_gb=$(df -BG --output=avail "${NXF_TASK_WORKDIR:-.}" | tail -1 | tr -dc '0-9')
if [ "$free_gb" -lt "$min_gb" ]; then
  echo "DISK GUARD: only ${free_gb} GB free (< ${min_gb} GB). Clear completed work directories before continuing." >&2
  exit 75
fi
