#!/usr/bin/env bash
# Periodic status snapshot for the long-running S1.5 / S2.1 jobs.
# Appends to deposits/jobs_status.log every 30 minutes; negligible CPU.
# Started with:  nohup bash monitor_jobs.sh > /dev/null 2>&1 &
cd "$(dirname "$0")" || exit 1

snap() {
  echo "=== $(date '+%Y-%m-%d %H:%M') ==============================================="
  for f in deposits/env/w5rig_k5.log deposits/env/w5rig_k8.log \
    deposits/ar2/ar1_r2000b.log deposits/ar2/effdf.log; do
    if [ -f "$f" ]; then
      printf '%-36s %8s bytes | ' "$f" "$(wc -c <"$f" | tr -d ' ')"
      tail -1 "$f" 2>/dev/null | head -c 200
      echo
    fi
  done
  printf 'python procs: %s | ' "$(powershell -NoProfile -Command \
    "(Get-CimInstance Win32_Process -Filter \"name = 'python.exe'\" | Measure-Object).Count" \
    2>/dev/null | tr -d '\r' | tail -1)"
  printf 'total CPU: %s\n' "$(powershell -NoProfile -Command \
    "'{0:N0}%' -f (Get-Counter '\Processor(_Total)\% Processor Time' -SampleInterval 1 -MaxSamples 1).CounterSamples.CookedValue" \
    2>/dev/null | tr -d '\r' | tail -1)"
  # heavy workers still burning (script name is absent from worker argv)
  powershell -NoProfile -Command \
    "Get-CimInstance Win32_Process -Filter \"name = 'python.exe'\" | ForEach-Object { \$pr=Get-Process -Id \$_.ProcessId -EA SilentlyContinue; if (\$pr.CPU -gt 60) { '  worker {0} cpu={1:N0}s' -f \$_.ProcessId, \$pr.CPU } }" \
    2>/dev/null | tr -d '\r'
}

for _ in $( # ~48 hours of coverage
  seq 1 96
); do
  snap >>deposits/jobs_status.log 2>&1
  sleep 1800
done
