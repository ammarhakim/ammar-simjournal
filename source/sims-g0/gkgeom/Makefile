# -*- makefile-gmake -*-

# Sample Makefile to use installed gkylzero library: copy and modify
# for your needs

# Default lapack include and libraries: we prefer linking to static library
LAPACK_INC = ${HOME}/gkylsoft/OpenBLAS/include
LAPACK_LIB_DIR = ${HOME}/gkylsoft/OpenBLAS/lib
LAPACK_LIB = -lopenblas

# SuperLU includes and librararies
SUPERLU_INC = ${HOME}/gkylsoft/superlu/include
ifeq ($(UNAME_S),Linux)
	SUPERLU_LIB_DIR = ${HOME}/gkylsoft/superlu/lib64
	SUPERLU_LIB = ${HOME}/gkylsoft/superlu/lib64/libsuperlu.a
else
	SUPERLU_LIB_DIR = ${HOME}/gkylsoft/superlu/lib
	SUPERLU_LIB = ${HOME}/gkylsoft/superlu/lib/libsuperlu.a
endif

PREFIX = ${HOME}/gkylsoft

# Include config.mak file (if it exists)
-include config.mak

CFLAGS = -O3 -g -ffast-math -I. -fPIC

G0_INC_DIR = ${PREFIX}/gkylzero/include
G0_LIB_DIR = ${PREFIX}/gkylzero/lib
G0_LIB = -lgkylzero
CUDA_LIB = 

ifeq ($(CC), nvcc)
	CUDA_LIB = -lcublas
endif

G0_LIBS = ${G0_LIB} ${CUDA_LIB} -lm -lpthread
G0_RPATH = -Wl,-rpath,${G0_LIB_DIR}

# determine OS we are running on
UNAME = $(shell uname)

EXTRA_LINK_FLAGS = 

# On OSX we should use Accelerate framework
ifeq ($(UNAME), Darwin)
	LAPACK_INC = . # dummy
	LAPACK_LIB = -framework Accelerate
	CFLAGS += -DGKYL_USING_FRAMEWORK_ACCELERATE
	EXTRA_LINK_FLAGS += -pagezero_size 10000 -image_base 100000000
endif

# On Linux we need to pass proper link flags to export all symbols
ifeq ($(UNAME), Linux)
	EXTRA_LINK_FLAGS += -Wl,-E
endif

# Read MPI paths and flags if needed 
USING_MPI =
MPI_INC_DIR = . # dummy
MPI_LIB_DIR = .
ifeq (${USE_MPI}, 1)
	USING_MPI = yes
	MPI_INC_DIR = ${MPI_INC}
	MPI_LIB_DIR = ${MPI_LIB}
	MPI_LIBS = -lmpi
	CFLAGS += -DGKYL_HAVE_MPI
endif

# Read LUA paths and flags if needed 
USING_LUA =
LUA_INC_DIR = . # dummy
LUA_LIB_DIR = .
ifeq (${CONF_USE_LUA}, 1)
	USING_LUA = yes
	LUA_INC_DIR = ${CONF_LUA_INC}
	LUA_LIB_DIR = ${CONF_LUA_LIB_DIR}
	LUA_LIBS = ${CONF_LUA_LIBS}
	CFLAGS += -DGKYL_HAVE_LUA
endif

INCLUDES = -I${G0_INC_DIR} -I${LAPACK_INC} -I${SUPERLU_INC} -I${MPI_INC_DIR} -I${LUA_INC_DIR}
LIB_DIRS = -L${LAPACK_LIB_DIR} -L${SUPERLU_LIB_DIR} -L${MPI_LIB_DIR} -L${LUA_LIB_DIR}
EXT_LIBS = ${LAPACK_LIB} ${SUPERLU_LIB} ${MPI_LIBS} 

all: geotest1 geotest2 geotest3 geotest4 rid-test ## build all targets

gkgeom.o: gkgeom.c gkgeom.h
	 ${CC} ${CFLAGS} ${INCLUDES} -c gkgeom.c -o gkgeom.o

dbl-exp.o: dbl-exp.c
	 ${CC} ${CFLAGS} ${INCLUDES} -c dbl-exp.c -o dbl-exp.o

ridders.o: ridders.c
	 ${CC} ${CFLAGS} ${INCLUDES} -c ridders.c -o ridders.o

geotest1: geotest1.c gkgeom.o dbl-exp.o ridders.o
	${CC} ${CFLAGS} ${INCLUDES} geotest1.c -o geotest1 gkgeom.o dbl-exp.o ridders.o -L${G0_LIB_DIR} ${G0_RPATH} ${G0_LIBS} ${LIB_DIRS} ${EXT_LIBS}

geotest2: geotest2.c gkgeom.o dbl-exp.o ridders.o
	${CC} ${CFLAGS} ${INCLUDES} geotest2.c -o geotest2 gkgeom.o dbl-exp.o ridders.o -L${G0_LIB_DIR} ${G0_RPATH} ${G0_LIBS} ${LIB_DIRS} ${EXT_LIBS}

geotest3: geotest3.c gkgeom.o dbl-exp.o ridders.o
	${CC} ${CFLAGS} ${INCLUDES} geotest3.c -o geotest3 gkgeom.o dbl-exp.o ridders.o -L${G0_LIB_DIR} ${G0_RPATH} ${G0_LIBS} ${LIB_DIRS} ${EXT_LIBS}

geotest4: geotest4.c gkgeom.o dbl-exp.o ridders.o
	${CC} ${CFLAGS} ${INCLUDES} geotest4.c -o geotest4 gkgeom.o dbl-exp.o ridders.o -L${G0_LIB_DIR} ${G0_RPATH} ${G0_LIBS} ${LIB_DIRS} ${EXT_LIBS}

rid-test: rid-test.c gkgeom.o dbl-exp.o ridders.o
	${CC} ${CFLAGS} ${INCLUDES} rid-test.c -o rid-test gkgeom.o dbl-exp.o ridders.o -L${G0_LIB_DIR} ${G0_RPATH} ${G0_LIBS} ${LIB_DIRS} ${EXT_LIBS}

.PHONY: clean
clean: ## Clean build outputs
	rm -rf gkgeom.o geotest1 geotest2 geotest3 geotest4

# From: https://www.client9.com/self-documenting-makefiles/
.PHONY: help
help: ## Show help
	@echo "Following targets are defined:"
	@awk -F ':|##' '/^[^\t].+?:.*?##/ {\
        printf "\033[36m%-30s\033[0m %s\n", $$1, $$NF \
        }' $(MAKEFILE_LIST)
