# Pipeline steps, in dependency order. Everything runs through `uv run`, which is what
# keeps four machines on the same interpreter and the same library versions.
#
# The ordering here is not decoration: feats.py must precede everything, and ablate.py
# must precede anything that reads results/preds/oof.json.

UV := uv run

.DEFAULT_GOAL := help
.PHONY: help setup hooks features baseline ablate score submit reweight verify verify-extra test doc clean-cache

help:  ## show this help
	@grep -hE '^[a-z-]+:.*?##' $(MAKEFILE_LIST) | sort | \
	 awk 'BEGIN{FS=":.*?## "}{printf "  \033[36m%-12s\033[0m %s\n", $$1, $$2}'

hooks:  ## arm the local pre-push guard on main (once per clone)
	@test -x .githooks/pre-push || { echo "no executable .githooks/pre-push"; exit 1; }
	git config core.hooksPath .githooks
	@echo "core.hooksPath = $$(git config --get core.hooksPath)"

setup: hooks  ## environment, submodule, data and features — run once per machine
	uv sync
	git submodule update --init
	bash data/fetch.sh
	$(UV) python src/feats.py
	@echo
	@echo "Готово. Проверь установку:  uv run pytest   (должно быть 6 passed)"

data/feats.npz: src/feats.py
	$(UV) python src/feats.py

features: data/feats.npz  ## build data/feats.npz and data/rows.csv (~45 s)

baseline: data/feats.npz  ## standalone baseline run (~10 min)
	$(UV) python src/base.py

ablate: data/feats.npz  ## regenerate results/preds/oof.json — SLOW (~1 h), see CONTRIBUTING.md
	$(UV) python src/ablate.py

score: data/feats.npz  ## metrics and paired bootstrap over the saved predictions (~2 min)
	$(UV) python src/score.py

submit: data/feats.npz  ## build both submission files and run the organisers' validators
	$(UV) python src/submit.py

reweight: data/feats.npz  ## score under a test-like label marginal (~3 min)
	$(UV) python src/reweight.py

trunk: data/feats.npz  ## both arms of the joint likelihood — SLOW (~1 h), writes results/preds/trunk_*.json
	uv run python src/trunk.py --mode twohead    --seeds 0,1,2,3 --lams 0,0.3,1.0,3.0
	uv run python src/trunk.py --mode calibrated --seeds 0,1,2,3 --lams 0,0.3,1.0,3.0

trunk-noise: data/feats.npz  ## the noise ladder on the screening channel — SLOW (~1 h)
	for e in 0.5 1 2 4; do \
	  uv run python src/trunk.py --mode twohead    --seeds 0,1,2,3 --lams 3.0 --noise $$e; \
	  uv run python src/trunk.py --mode calibrated --seeds 0,1,2,3 --lams 3.0 --noise $$e; \
	done

trunk-score: data/feats.npz  ## read the saved trunk predictions: lambda response, dose curve, noise curve (~5 min)
	uv run python src/trunkscore.py
	uv run python src/trunkdose.py
	uv run python src/trunknoise.py

test:  ## golden-value guard on the cross-validation split
	$(UV) pytest

verify: data/feats.npz  ## the quick verification scripts (skips f3, f12: ~70 min combined)
	@for f in verify/f1_formulas.py verify/f2_splitnorm.py verify/f6_data.py \
	          verify/f8_mnar.py verify/f9_mccthr.py verify/f10_calib.py \
	          verify/g1_calib.py verify/g2_factor.py; do \
	  echo "=== $$f ==="; $(UV) python $$f || exit 1; \
	done

verify-extra: data/feats.npz  ## the h* and k* verification scripts (~25 min), logs into results/logs/
	@mkdir -p results/logs
	@for f in h1_geometry h2_tdi_alerts h3_alerts_delta k1_shrink k3_center k4_enrich \
	          k5_shift k6_shift1d k7_2d6shift k8_kernel k9_shape; do \
	  echo "=== $$f ==="; $(UV) python verify/$$f.py > results/logs/$$f.log 2>&1 || exit 1; \
	done
	@echo "logs in results/logs/"

doc:  ## rebuild docs/CYP — модель и данные.pdf (needs XeLaTeX + ParaType)
	bash docs/build.sh

clean-cache:  ## drop the clustering cache; it rebuilds itself
	rm -f data/clusters_*.npz data/folds_*.npy data/folds_*.npz
