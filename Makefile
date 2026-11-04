PYTHON ?= python3
.PHONY: all check figures tables
all: check tables figures
check:
	$(PYTHON) 04_figures_tables/check.py
	$(PYTHON) 02_benchmark/verify_archive.py
figures:
	$(PYTHON) 04_figures_tables/make_architecture.py
	$(PYTHON) 04_figures_tables/make_figure1.py
	$(PYTHON) 04_figures_tables/make_agreement_fitting_figures.py
	$(PYTHON) 04_figures_tables/make_supplement_figures.py
	$(PYTHON) 04_figures_tables/make_speedup_figure.py
	$(PYTHON) 04_figures_tables/make_memory_figure.py
	$(PYTHON) 04_figures_tables/make_report_figure.py
	$(PYTHON) 04_figures_tables/make_report_figure_compact.py
tables:
	$(PYTHON) 04_figures_tables/tables.py
