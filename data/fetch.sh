#!/usr/bin/env bash
# Данные соревнования с Hugging Face.
set -euo pipefail
cd "$(dirname "$0")"

REPO="openadmet/cyp-challenge-train-test"
FILES=(
  "cyp-challenge-TRAIN_inhibition.csv"
  "cyp-challenge-TRAIN_TDI.csv"
  "cyp-challenge-TRAIN_Emax.csv"
  "cyp-challenge-single-concentration-TRAIN.csv"
  "cyp-challenge-TEST-BLINDED.csv"
)

for f in "${FILES[@]}"; do
  if [ -f "$f" ]; then echo "уже есть: $f"; continue; fi
  echo "качаю: $f"
  curl -fSL -o "$f" "https://huggingface.co/datasets/${REPO}/resolve/main/${f}"
done

echo
echo "готово, строк по файлам:"
for f in "${FILES[@]}"; do printf "  %-48s %s\n" "$f" "$(($(wc -l < "$f") - 1))"; done
