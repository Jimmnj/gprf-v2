#!/usr/bin/env bash
# fetch_diskmass.sh -- download the DiskMass Survey arXiv sources whose LaTeX tables dms_fb_20261105.py parses
# (Martinsson et al. 2013a, DMS VI, A&A 557, A130, arXiv:1307.8130; 2013b, DMS VII, A&A 557, A131, arXiv:1308.0336).
# Not redistributed here. Run from the repository root: bash tools/fetch_diskmass.sh
set -eu
mkdir -p data/diskmass && cd data/diskmass
for id in 1307.8130 1308.0336; do
  [ -d "arxiv_${id}_src" ] && { echo "have $id"; continue; }
  curl -L --fail --retry 2 --connect-timeout 20 --max-time 180 -A "gprf-v2-fetch" -o "arxiv_${id}_src.tgz" "https://arxiv.org/e-print/$id"
  mkdir -p "arxiv_${id}_src" && tar -xzf "arxiv_${id}_src.tgz" -C "arxiv_${id}_src" && echo "unpacked $id"
done
ls arxiv_1308.0336_src/tables arxiv_1307.8130_src/tables
