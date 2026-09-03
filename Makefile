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

.PHONY: all clean no-omp lego-tools afc dislo test check test-large

all: lego-tools afc dislo

lego-tools:
	$(MAKE) -C lego-tools all

afc:
	$(MAKE) -C afc all

dislo: lego-tools
	$(MAKE) -C dislo all

no-omp:
	$(MAKE) -C lego-tools no-omp
	$(MAKE) -C afc OPENMP=no all
	$(MAKE) -C dislo no-omp

# Run the end-to-end smoke test. Builds first so tests see up-to-date binaries.
# `check` is a conventional alias.
test check: all
	@bash tests/smoke.sh

# Regression test for the large-file (>4 GiB) I/O path: rebuilds lego-tools
# twice with different GZ_CHUNK_MAX values and compares byte-identical
# output.  Slower than `make test` because it does two full rebuilds.
test-large: all
	@bash tests/large-file.sh

clean:
	$(MAKE) -C lego-tools clean
	$(MAKE) -C afc clean
	$(MAKE) -C dislo clean
