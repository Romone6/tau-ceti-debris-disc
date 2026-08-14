PYTHON ?= python3
export PYTHONPATH := $(CURDIR)/src

verify:
	$(PYTHON) scripts/verify_release.py

reproduce:
	$(PYTHON) scripts/reproduce_collisional.py

test:
	$(PYTHON) -m pytest -q
