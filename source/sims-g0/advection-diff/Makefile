# -*- makefile-gmake -*-

ARCH_FLAGS ?= -march=native
PREFIX ?= ${HOME}/gkylsoft
INSTALL_PREFIX ?= ${PREFIX}

# determine OS we are running on
UNAME = $(shell uname)

# Include config.mak file (if it exists) to overide defaults above
-include $(PREFIX)/gminus/share/config.mak

CFLAGS = -O3 -g -ffast-math 

G0_INC_DIR = ${PREFIX}/gminus/include
G0_LIB_DIR = ${PREFIX}/gminus/lib
G0_LIB = -lg0core -lm 

all: sim-1d-const-adiff sim-2d-rotflow sim-vlasov-pot

adiff.o: adiff.c adiff.h
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) -c adiff.c

sim-1d-const-adiff: sim-1d-const-adiff.c adiff.o
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) sim-1d-const-adiff.c adiff.o -L $(G0_LIB_DIR) $(G0_LIB) -o sim-1d-const-adiff

sim-2d-rotflow: sim-2d-rotflow.c adiff.o
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) sim-2d-rotflow.c adiff.o -L $(G0_LIB_DIR) $(G0_LIB) -o sim-2d-rotflow

sim-vlasov-pot: sim-vlasov-pot.c adiff.o
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) sim-vlasov-pot.c adiff.o -L $(G0_LIB_DIR) $(G0_LIB) -o sim-vlasov-pot

.PHONY: clean
clean:
	rm -f adiff.o sim-1d-const-adiff sim-2d-rotflow sim-vlasov-pot
