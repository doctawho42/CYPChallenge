#!/usr/bin/env bash
# Build the main document. Needs XeLaTeX and the ParaType fonts:
#   Debian/Ubuntu: apt install texlive-xetex texlive-lang-cyrillic fonts-paratype
#   macOS:         brew install --cask mactex font-pt-serif font-pt-sans
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$(pwd)

if [ ! -f "$ROOT/data/cyp-challenge-TRAIN_inhibition.csv" ]; then
  echo "no data: run  bash data/fetch.sh" >&2; exit 1
fi
if [ ! -f "$ROOT/data/rows.csv" ]; then
  echo "no data/rows.csv: run  uv run python src/feats.py" >&2; exit 1
fi
if [ ! -f "$ROOT/results/nn_seed0.npy" ]; then
  echo "no results/nn_seed0.npy: run  uv run python verify/f12_cvhard.py  (~10 min)" >&2; exit 1
fi

# Figures are computed from the data, so they come first. figs.py writes straight into
# docs/tex/fig/, which is where the .tex files include them from.
uv run python docs/tex/figs.py

cd "$ROOT/docs/tex"
# Без этого XeLaTeX штампует в PDF время сборки, и закоммиченный артефакт меняется при
# каждой пересборке, ничего не меняющей по существу. Файл бинарный и не мержится, а машин
# четыре, так что каждый такой штамп --- заготовка конфликта на ровном месте. Дата взята
# фиксированной, а не текущей: смысл в том, чтобы одинаковый источник давал одинаковый файл.
export SOURCE_DATE_EPOCH=1700000000
export FORCE_SOURCE_DATE=1
for i in 1 2 3; do xelatex -interaction=nonstopmode main.tex >/dev/null; done
mv main.pdf "../CYP — модель и данные.pdf"
echo "done: docs/CYP — модель и данные.pdf"
