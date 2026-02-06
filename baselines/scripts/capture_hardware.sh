#!/usr/bin/env bash
set -euo pipefail
mkdir -p baselines
out="baselines/hardware.txt"
{
  echo "== date =="; date || true
  echo
  echo "== uname =="; uname -a || true
  echo
  echo "== os (linux) =="; (command -v lsb_release >/dev/null 2>&1 && lsb_release -a) || true
  echo
  echo "== os (mac) =="; (command -v sw_vers >/dev/null 2>&1 && sw_vers) || true
  echo
  echo "== cpu =="; (command -v lscpu >/dev/null 2>&1 && lscpu) || (command -v sysctl >/dev/null 2>&1 && sysctl -n machdep.cpu.brand_string) || true
  echo
  echo "== mem =="; (command -v free >/dev/null 2>&1 && free -h) || (command -v sysctl >/dev/null 2>&1 && sysctl -n hw.memsize) || true
  echo
  echo "== gpu (nvidia) =="; (command -v nvidia-smi >/dev/null 2>&1 && nvidia-smi) || true
} > "$out"
echo "Wrote $out"
