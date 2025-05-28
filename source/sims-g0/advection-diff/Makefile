# -*- makefile-gmake -*-

ARCH_FLAGS ?= -march=native
PREFIX ?= ${HOME}/gkylsoft

# Include config.mak file (if it exists) to overide defaults above
-include $(PREFIX)/gminus/share/config.mak

CFLAGS = -O3 -g -ffast-math 

G0_INC_DIR = ${PREFIX}/gminus/include
G0_LIB_DIR = ${PREFIX}/gminus/lib
G0_LIB = -lg0core -lm 

all: sim-1d-sin sim-1d-gauss sim-2d-rotflow sim-vlasov-cos-pot sim-vlasov-shm-pot

adiff.o: adiff.c adiff.h
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) -c adiff.c

sim-1d-sin: sim-1d-sin.c adiff.o
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) sim-1d-sin.c adiff.o -L $(G0_LIB_DIR) $(G0_LIB) -o sim-1d-sin

sim-1d-gauss: sim-1d-gauss.c adiff.o
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) sim-1d-gauss.c adiff.o -L $(G0_LIB_DIR) $(G0_LIB) -o sim-1d-gauss

sim-2d-rotflow: sim-2d-rotflow.c adiff.o
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) sim-2d-rotflow.c adiff.o -L $(G0_LIB_DIR) $(G0_LIB) -o sim-2d-rotflow

sim-vlasov-cos-pot: sim-vlasov-cos-pot.c adiff.o
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) sim-vlasov-cos-pot.c adiff.o -L $(G0_LIB_DIR) $(G0_LIB) -o sim-vlasov-cos-pot

sim-vlasov-shm-pot: sim-vlasov-shm-pot.c adiff.o
	$(CC) $(CFLAGS) $(ARCH_FLAGS) -I. -I $(G0_INC_DIR) sim-vlasov-shm-pot.c adiff.o -L $(G0_LIB_DIR) $(G0_LIB) -o sim-vlasov-shm-pot

.PHONY: clean
clean:
	rm -f adiff.o sim-1d-gauss sim-1d-sin sim-2d-rotflow sim-vlasov-cos-pot sim-vlasov-shm-pot
