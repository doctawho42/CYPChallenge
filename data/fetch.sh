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

# Внешние метки: обучающий набор базовой модели CheMeleon, опубликованной организаторами.
# 8068 соединений, Apache-2.0, pIC50 по тем же четырём ферментам, курировано из ChEMBL.
# Почему это записано здесь, а не подразумевается: пункты 60-65, 77, 144 и 290 стоят на этом
# наборе, и до 23 сентября НИГДЕ в дереве не было сказано, откуда он берётся -- ни адреса, ни
# скрипта, ни суммы (пункт 318). Файл, который нельзя переполучить, -- это результат, который
# нельзя воспроизвести.
#
# Берём репозиторий -chemeleon-baseline, а не -chemeleon-v1: тела файлов у них ПОБАЙТОВО
# одинаковы, различаются только заголовки регистром, и src/trunkext.py с src/ablext.py читают
# OPENADMET_LOGAC50_{нижний регистр}, то есть вариант baseline.
EXT_REPO="openadmet/cyp1a2-cyp2d6-cyp3a4-cyp3c9-chemeleon-baseline"
declare -a EXT_SRC=("X_train.csv" "y_train.csv")
declare -a EXT_DST=("chemeleon_X_train.csv" "chemeleon_y_train.csv")
declare -a EXT_SHA=(
  "208d3662cdf3221c9ef7e4adc138fd48e1ae03ef9edb3e984f8344c8c26d35dd"
  "1ab43020a8064c4dac065afa60f65d021692f1716a4aacec1b6808c44fb37d9e"
)
for i in "${!EXT_SRC[@]}"; do
  d="${EXT_DST[$i]}"
  if [ ! -f "$d" ]; then
    echo "качаю: $d"
    curl -fSL -o "$d" \
      "https://huggingface.co/${EXT_REPO}/resolve/main/anvil_training/data/${EXT_SRC[$i]}"
  else
    echo "уже есть: $d"
  fi
  # Сумма сверяется ВСЕГДА, а не только после скачивания: файл на диске мог приехать откуда
  # угодно, и провенанс -- это про то, что лежит сейчас, а не про то, что однажды качали.
  got=$(shasum -a 256 "$d" | cut -d' ' -f1)
  if [ "$got" != "${EXT_SHA[$i]}" ]; then
    echo "ОСТАНОВ: $d не та, что записана здесь." >&2
    echo "  ожидалось ${EXT_SHA[$i]}" >&2
    echo "  получено  $got" >&2
    echo "  Внешние метки входят в обучение, и подменённый файл сдвинул бы каждое число." >&2
    exit 1
  fi
done

echo
echo "готово, строк по файлам:"
for f in "${FILES[@]}" "${EXT_DST[@]}"; do printf "  %-48s %s\n" "$f" "$(($(wc -l < "$f") - 1))"; done
