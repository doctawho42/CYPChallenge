#!/usr/bin/env bash
# Standalone build of the biologist-facing document. Needs XeLaTeX and ParaType fonts,
# the same toolchain as docs/build.sh, but no data and no figures.
set -euo pipefail
cd "$(dirname "$0")"
export SOURCE_DATE_EPOCH=1725840000
xelatex -interaction=nonstopmode -halt-on-error bio.tex >/dev/null
xelatex -interaction=nonstopmode -halt-on-error bio.tex >/dev/null
echo "docs/bio/bio.pdf"
