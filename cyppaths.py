"""Пути внутри репозитория. Скрипты подключают это вместо жёстких путей.

DATA     — данные соревнования, каталог data/
RESULTS  — логи и сохранённые предсказания, каталог results/
TUTORIAL — репозиторий оргкомитета с официальной реализацией метрики.
           Ищется рядом с этим репозиторием либо по переменной окружения CYP_TUTORIAL:
             git clone https://github.com/OpenADMET/CYP-Challenge-Tutorial
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
PREDS = RESULTS / "preds"
LOGS = RESULTS / "logs"
DOCS = ROOT / "docs"

# многие скрипты писались со строковым префиксом, поэтому оставляем и такую форму
D = str(DATA) + os.sep
RES = str(RESULTS) + os.sep


def tutorial() -> Path:
    """Каталог CYP-Challenge-Tutorial; кладёт его в sys.path и возвращает."""
    env = os.environ.get("CYP_TUTORIAL")
    candidates = [Path(env)] if env else []
    candidates += [ROOT.parent / "CYP-Challenge-Tutorial", ROOT / "CYP-Challenge-Tutorial"]
    for c in candidates:
        if (c / "evaluation" / "custom_scoring_functions.py").exists():
            if str(c) not in sys.path:
                sys.path.insert(0, str(c))
            return c
    raise FileNotFoundError(
        "не найден CYP-Challenge-Tutorial. Клонируйте его рядом с репозиторием:\n"
        "  git clone https://github.com/OpenADMET/CYP-Challenge-Tutorial\n"
        "или укажите путь в переменной окружения CYP_TUTORIAL"
    )
