
# Sample Makefile to use installed gkylzero library: copy and modify
# for your needs

CFLAGS = -O3 -g -I.
PREFIX = ${HOME}/gkylsoft

G0_INC_DIR = ${PREFIX}/gkylzero/include
G0_LIB_DIR = ${PREFIX}/gkylzero/lib
G0_LIB_FLAGS = -lgkylzero -lm

all: kh-euler

kh-euler: kh-euler.c
	 ${CC} ${CFLAGS} -I${G0_INC_DIR} kh-euler.c -o kh-euler -L${G0_LIB_DIR} ${G0_LIB_FLAGS} 

clean:
	rm -rf kh-euler kh-euler.dSYM

