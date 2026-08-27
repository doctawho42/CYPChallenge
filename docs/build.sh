#!/usr/bin/env bash
# Сборка основного документа. Нужен XeLaTeX и шрифты ParaType.
#   Debian/Ubuntu: apt install texlive-xetex texlive-lang-cyrillic fonts-paratype
#   macOS:         brew install --cask mactex && brew install --cask font-pt-serif font-pt-sans
set -euo pipefail
cd "$(dirname "$0")/tex"

# графики строятся из данных, поэтому сначала данные
if [ ! -f ../../data/cyp-challenge-TRAIN_inhibition.csv ]; then
  echo "нет данных, запустите data/fetch.sh" >&2; exit 1
fi
python3 figs.py

for i in 1 2 3; do xelatex -interaction=nonstopmode main.tex >/dev/null; done
mv main.pdf "../CYP — модель и данные.pdf"
echo "готово: docs/CYP — модель и данные.pdf"
