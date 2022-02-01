EIGEN_INCLUDES=/Users/ahakim/gkylsoft/eigen/include/eigen3

all: dg-int

dg-int: dg-int.cxx
	g++ -O3 dg-int.cxx -I${EIGEN_INCLUDES} -o dg-int

clean:
	rm -rf dg-int
