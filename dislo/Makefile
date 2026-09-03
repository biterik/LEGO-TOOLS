# dislo/Makefile — builds lego-dislo
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

ATOMIO = ../lego-tools/atomio.o

all: lego-dislo

$(ATOMIO): ../lego-tools/atomio.c ../lego-tools/atomio.h
	$(MAKE) -C ../lego-tools atomio.o

lego-dislo: main.c card.c sextic.c findroots.c \
            card.h sextic.h ../lego-tools/atomio.h $(ATOMIO)
	$(CC) $(CFLAGS) main.c card.c sextic.c findroots.c $(ATOMIO) \
	    -o lego-dislo $(LDFLAGS) $(LDLIBS)

# Build without OpenMP
no-omp:
	$(MAKE) OMP_FLAG= OMP_INC= OMP_LIB=

clean:
	rm -f lego-dislo

.PHONY: all clean no-omp
