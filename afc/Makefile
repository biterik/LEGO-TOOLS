# afc/Makefile — Atomic Format Converter
#
# Author:  Erik Bitzek <e.bitzek@mpi-susmat.de>
#          Max-Planck-Institut fuer Nachhaltige Materialien,
#          Duesseldorf, Germany
# Funding: NFDI-MatWerk

CC       ?= cc
CFLAGS   ?= -O2 -Wall -Wextra -std=c11
LDFLAGS  ?=
LDLIBS    = -lz -lm

# OpenMP — disable with: make OPENMP=no
OPENMP   ?= yes
ifeq ($(OPENMP),yes)
  UNAME_S := $(shell uname -s)
  ifeq ($(UNAME_S),Darwin)
    LIBOMP_PREFIX := $(shell brew --prefix libomp 2>/dev/null)
    ifneq ($(LIBOMP_PREFIX),)
      CFLAGS  += -Xpreprocessor -fopenmp -I$(LIBOMP_PREFIX)/include
      LDFLAGS += -L$(LIBOMP_PREFIX)/lib
      LDLIBS  += -lomp
    endif
  else
    CFLAGS  += -fopenmp
    LDFLAGS += -fopenmp
  endif
endif

TARGET = afc
SRCS   = main.c converter.c
HDRS   = converter.h elements.h

all: $(TARGET)

$(TARGET): $(SRCS) $(HDRS)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(SRCS) $(LDLIBS)

clean:
	rm -f $(TARGET)

.PHONY: all clean
