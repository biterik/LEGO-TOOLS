# LEGO-TOOLS — top-level Makefile
#
# Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
#          Max-Planck-Institut fuer Nachhaltige Materialien,
#          Duesseldorf, Germany
# Funding: NFDI-MatWerk
#
# Builds both the lego-tools suite and the Atomic Format Converter (afc).
# Run with:
#   make          — build everything
#   make no-omp   — build without OpenMP (e.g. if libomp is missing)
#   make clean    — remove all build artefacts

.PHONY: all clean no-omp lego-tools afc test check

all: lego-tools afc

lego-tools:
	$(MAKE) -C lego-tools all

afc:
	$(MAKE) -C afc all

no-omp:
	$(MAKE) -C lego-tools no-omp
	$(MAKE) -C afc OPENMP=no all

# Run the end-to-end smoke test. Builds first so tests see up-to-date binaries.
# `check` is a conventional alias.
test check: all
	@bash tests/smoke.sh

clean:
	$(MAKE) -C lego-tools clean
	$(MAKE) -C afc clean
