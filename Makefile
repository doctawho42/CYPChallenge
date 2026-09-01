# Pipeline steps, in dependency order. Everything runs through `uv run`, which is what
# keeps four machines on the same interpreter and the same library versions.
#
# The ordering here is not decoration: feats.py must precede everything, and ablate.py
# must precede anything that reads results/preds/oof.json.

UV := uv run

.DEFAULT_GOAL := help
.PHONY: help setup hooks features baseline ablate score submit reweight verify \
        verify-extra verify-delta verify-regime test doc clean-cache

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

test:  ## golden-value guard on the split, plus the fifth ensemble member's guards
	$(UV) pytest

verify: data/feats.npz  ## the quick verification scripts (skips f3, f12: ~70 min combined)
	@for f in verify/f1_formulas.py verify/f2_splitnorm.py verify/f6_data.py \
	          verify/f8_mnar.py verify/f9_mccthr.py verify/f10_calib.py \
	          verify/g1_calib.py verify/g2_factor.py; do \
	  echo "=== $$f ==="; $(UV) python $$f || exit 1; \
	done

verify-extra: data/feats.npz  ## h* and k1-k10: rescaling and how far the test sits (~25 min)
	@mkdir -p results/logs
	@for f in h1_geometry h2_tdi_alerts h3_alerts_delta k1_shrink k3_center k4_enrich \
	          k5_shift k6_shift1d k7_2d6shift k8_kernel k9_shape k10_strat2d6; do \
	  echo "=== $$f ==="; $(UV) python verify/$$f.py > results/logs/$$f.log 2>&1 || exit 1; \
	done
	@echo "logs in results/logs/"

verify-delta: data/feats.npz  ## k11-k19: external data and the test label shift (~55 min)
	@mkdir -p results/logs
	@for f in k11_exttransfer k12_extneighbors k13_channels k14_design k15_pooldelta \
	          k16_modelspread k17_ensdelta k18_nbspace k19_ens3delta; do \
	  echo "=== $$f ==="; $(UV) python verify/$$f.py > results/logs/$$f.log 2>&1 || exit 1; \
	done
	@echo "logs in results/logs/"

verify-regime: data/feats.npz  ## k20-k27: is our regime the test's, and what it costs (~40 min)
	@mkdir -p results/logs
	@for f in k20_strat k21_borda k22_layerboot k23_tilt k24_visible k25_reweight \
	          k26_screen k27_trunkens; do \
	  echo "=== $$f ==="; $(UV) python verify/$$f.py > results/logs/$$f.log 2>&1 || exit 1; \
	done
	@echo "logs in results/logs/"

verify-bounds: data/feats.npz  ## k28-k35: what is provably out of reach, and the anchor split (~50 min)
	@mkdir -p results/logs
	@for f in k28_ess k29_positivity k30_oracle k31_campaign k32_anchor k33_lbident \
	          k34_series k35_plate; do \
	  echo "=== $$f ==="; $(UV) python verify/$$f.py > results/logs/$$f.log 2>&1 || exit 1; \
	done
	@echo "logs in results/logs/"

verify-ceiling: data/feats.npz  ## k36-k44: the gap to the screen, the criterion, the learner (~60 min)
	@mkdir -p results/logs
	@for f in k36_ceiling k37_gap k38_trunc k39_splits k40_topk k41_earlystop \
	          k42_visiblerank k43_lbpredict k44_bits; do \
	  echo "=== $$f ==="; $(UV) python verify/$$f.py > results/logs/$$f.log 2>&1 || exit 1; \
	done
	@echo "logs in results/logs/"

ncgc:  ## fetch, merge and featurise the NCGC panel (~25 min, network on first run)
	$(UV) python src/fetch1851.py
	$(UV) python src/ncgcmerge.py
	$(UV) python src/ncgcfeats.py

aux: data/feats.npz  ## the screening table as a training target — SLOW (~4 h)
	$(UV) python src/ablaux.py --seeds 0

ncgc-ablate: data/feats.npz  ## the NCGC panel as extra training rows — SLOW (~10 h)
	$(UV) python src/ablncgc.py --seeds 0

deadpair: data/feats.npz  ## do the dead zone and the pairwise loss add, or overlap (~80 min)
	$(UV) python src/abldeadpair.py --seeds 0,1,2,3

doc:  ## rebuild docs/CYP — модель и данные.pdf (needs XeLaTeX + ParaType)
	bash docs/build.sh

clean-cache:  ## drop the clustering cache; it rebuilds itself
	rm -f data/clusters_*.npz data/folds_*.npy data/folds_*.npz
