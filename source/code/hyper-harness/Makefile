# -*- makefile-gmake -*-

# Makefile for C input files located outside of the gkeyll/ repository.

CC = clang
CFLAGS ?= -O3 -g -ffast-math -fPIC -MMD -MP
.PHONY: all clean help moments moments-unit moments-regression

all: ctest_eqn_euler ctest_eqn_maxwell ctest_eqn_burgers ctest_eqn_isoeuler ctest_eqn_shallow_water test_app

eqn_headers:
	maxima -b burgers.mac
	maxima -b euler.mac
	maxima -b isoeuler.mac
	maxima -b maxwell.mac
	maxima -b shallow-water.mac

mug.o: mug.c
	${CC} ${CFLAGS} -c mug.c

test_app: test_app.c app.h app.c mug.o
	${CC} ${CFLAGS} -I.  mug.o app.c $< -o $@ -lm -lpthread

ctest_eqn_euler: ctest_eqn_euler.c eqn_euler.h mug.o
	${CC} ${CFLAGS} -I. mug.o  $< -o $@ -lm -lpthread

ctest_eqn_maxwell: ctest_eqn_maxwell.c eqn_maxwell.h mug.o
	${CC} ${CFLAGS} -I. mug.o  $< -o $@ -lm -lpthread

ctest_eqn_burgers: ctest_eqn_burgers.c eqn_burgers.h mug.o
	${CC} ${CFLAGS} -I. mug.o $< -o $@ -lm -lpthread

ctest_eqn_isoeuler: ctest_eqn_isoeuler.c eqn_isoeuler.h mug.o
	${CC} ${CFLAGS} -I. mug.o  $< -o $@ -lm -lpthread

ctest_eqn_shallow_water: ctest_eqn_shallow_water.c eqn_shallow_water.h mug.o
	${CC} ${CFLAGS} -I. mug.o $< -o $@ -lm -lpthread

clean:
	rm -rf test_app ctest_eqn_euler ctest_eqn_maxwell ctest_eqn_isoeuler ctest_eqn_shallow_water ctest_eqn_burgers
