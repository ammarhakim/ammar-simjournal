
# Sample Makefile to use installed gkylzero library: copy and modify
# for your needs

CFLAGS = -O3 -g -I.
PREFIX = ${HOME}/gkylsoft

G0_INC_DIR = ${PREFIX}/gkylzero/include
G0_LIB_DIR = ${PREFIX}/gkylzero/lib
G0_LIB_FLAGS = -lm -lgkylzero

all: kh-euler

kh-euler: kh-euler.c
	 ${CC} ${CFLAGS} -I${G0_INC_DIR} -L${G0_LIB_DIR} ${G0_LIB_FLAGS} kh-euler.c -o kh-euler

clean:
	rm -rf kh-euler kh-euler.dSYM

