PYTHON ?= python3
OUTPUT ?= results/SUPP_TABLES_BRCA12_reproduced.xlsx
FIGURE_PREFIX ?= figures/reproduced/supp_fig2
EVE_ARGS ?=

.PHONY: install data reproduce verify test

install:
	$(PYTHON) -m pip install -r requirements-lock.txt

data:
	$(PYTHON) scripts/build_eve_artifacts.py $(EVE_ARGS)

reproduce: data
	$(PYTHON) main.py --output-workbook $(OUTPUT) --figure-prefix $(FIGURE_PREFIX)

verify:
	PYTHONPATH=src $(PYTHON) scripts/verify_reproduction.py $(OUTPUT) --figure-prefix $(FIGURE_PREFIX)

test:
	PYTHONPATH=src $(PYTHON) -m pytest
