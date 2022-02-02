
# Sample Makefile to use installed gkylzero library: copy and modify
# for your needs

CFLAGS = -O3 -g
PREFIX = ${HOME}/gkylsoft
ARCH_FLAGS = -march=native -ffast-math

G0_INC_DIR = ${PREFIX}/gkylzero/include
G0_SHARE_DIR = ${PREFIX}/gkylzero/share
G0_LIB_DIR = ${PREFIX}/gkylzero/lib
G0_LIB_FLAGS = -lgkylzero -lm

LAPACK_INC = ${HOME}/gkylsoft/OpenBLAS/include
LAPACK_LIB = ${HOME}/gkylsoft/OpenBLAS/lib/libopenblas.a

INCLUDES = -I${G0_INC_DIR} -I${LAPACK_INC} -I${G0_SHARE_DIR}

NVCC = 
USING_NVCC =
NVCC_FLAGS = 
CUDA_LIBS = 
ifeq ($(CC), nvcc)
       CFLAGS = -O3 -g --forward-unknown-to-host-compiler
       USING_NVCC = yes
       NVCC_FLAGS = -x cu -dc -arch=sm_70 --compiler-options="-fPIC" 
       LDFLAGS += -arch=sm_70
       CUDA_LIBS = -lcublas
endif

all: shock

shock: shock.c
	 ${CC} ${CFLAGS} $(INCLUDES) shock.c -o shock -L${G0_LIB_DIR} ${G0_LIB_FLAGS} ${LAPACK_LIB} ${CUDA_LIBS} -lpthread

clean:
	rm -rf shock

