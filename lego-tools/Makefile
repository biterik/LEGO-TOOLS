# lego-tools/Makefile
#
# Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
#          Max-Planck-Institut fuer Nachhaltige Materialien,
#          Duesseldorf, Germany
# Funding: NFDI-MatWerk

CC      ?= cc
CFLAGS  ?= -O3 -Wall -Wextra -std=c11 -g
LDFLAGS ?=
LDLIBS  ?= -lz -lm

# OpenMP: auto-detect on macOS (Apple clang needs libomp from Homebrew)
UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Darwin)
  OMP_INC  ?= -I/opt/homebrew/opt/libomp/include
  OMP_LIB  ?= -L/opt/homebrew/opt/libomp/lib -lomp
  OMP_FLAG := -Xpreprocessor -fopenmp
else
  OMP_INC  ?=
  OMP_LIB  ?= -lgomp
  OMP_FLAG := -fopenmp
endif

CFLAGS += $(OMP_FLAG) $(OMP_INC)
LDLIBS += $(OMP_LIB)

BINS = lego-cut lego-shift lego-analyze lego-remove-per-atom \
       lego-change-box lego-nearest-atoms lego-pbc-wrap lego-make-box

all: $(BINS)

atomio.o: atomio.c atomio.h
	$(CC) $(CFLAGS) -c atomio.c -o atomio.o

lego-cut: cut.c atomio.o atomio.h
	$(CC) $(CFLAGS) cut.c atomio.o -o lego-cut $(LDFLAGS) $(LDLIBS)

lego-shift: shift.c atomio.o atomio.h
	$(CC) $(CFLAGS) shift.c atomio.o -o lego-shift $(LDFLAGS) $(LDLIBS)

lego-analyze: analyze.c atomio.o atomio.h
	$(CC) $(CFLAGS) analyze.c atomio.o -o lego-analyze $(LDFLAGS) $(LDLIBS)

lego-remove-per-atom: remove_per_atom.c atomio.o atomio.h
	$(CC) $(CFLAGS) remove_per_atom.c atomio.o -o lego-remove-per-atom $(LDFLAGS) $(LDLIBS)

lego-change-box: change_box.c
	$(CC) $(CFLAGS) change_box.c -o lego-change-box $(LDFLAGS) $(LDLIBS)

lego-nearest-atoms: nearest_atoms.c atomio.o atomio.h
	$(CC) $(CFLAGS) nearest_atoms.c atomio.o -o lego-nearest-atoms $(LDFLAGS) $(LDLIBS)

lego-pbc-wrap: pbc_wrap.c atomio.o atomio.h
	$(CC) $(CFLAGS) pbc_wrap.c atomio.o -o lego-pbc-wrap $(LDFLAGS) $(LDLIBS)

lego-make-box: make_box.c atomio.o atomio.h
	$(CC) $(CFLAGS) make_box.c atomio.o -o lego-make-box $(LDFLAGS) $(LDLIBS)

# Build without OpenMP
no-omp:
	$(MAKE) OMP_FLAG= OMP_INC= OMP_LIB=

clean:
	rm -f $(BINS) atomio.o

.PHONY: all clean no-omp
