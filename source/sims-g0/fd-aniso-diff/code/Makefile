# -*- makefile -*-

# Sample Makefile to use installed gkylzero library: copy and modify
# for your needs

CFLAGS = -O3 -g -I. 
PREFIX = ${HOME}/gkylsoft

G0_INC_DIR = ${PREFIX}/gkylzero/include
G0_LIB_DIR = ${PREFIX}/gkylzero/lib
G0_LIBS = -lgkylzero -lm -lpthread

# Default lapack include and libraries: we prefer linking to static library
LAPACK_INC = ${HOME}/gkylsoft/OpenBLAS/include
LAPACK_LIB = ${HOME}/gkylsoft/OpenBLAS/lib/libopenblas.a

# determine OS we are running on
UNAME = $(shell uname)

# On OSX we should use Accelerate framework
ifeq ($(UNAME), Darwin)
	LAPACK_INC = . # dummy
	LAPACK_LIB = -framework Accelerate
	CFLAGS += -DGKYL_USING_FRAMEWORK_ACCELERATE
endif

INCLUDES = -I${G0_INC_DIR} -I${LAPACK_INC}

all: fd-aniso-diff

fd-aniso-diff: fd-aniso-diff.c
	 ${CC} ${CFLAGS} ${INCLUDES} fd-aniso-diff.c -o fd-aniso-diff -L${G0_LIB_DIR} ${G0_LIBS} ${LAPACK_LIB}

clean:
	rm -rf fd-aniso-diff fd-aniso-diff.dSYM

