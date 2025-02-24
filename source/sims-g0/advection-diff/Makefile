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

all: sim-1d-const-adiff 

adiff.o: adiff.c
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) -c adiff.c

sim-1d-const-adiff: sim-1d-const-adiff.c adiff.o
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) sim-1d-const-adiff.c -L $(G0_LIB_DIR) $(G0_LIB) -o sim-1d-const-adiff

.PHONY: clean
clean:
	rm -f sim-1d-const-adiff
