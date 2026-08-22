all: ctest_ga3

ga3_generated.h: ga3.mac
	maxima -b ga3.mac

ctest_ga3: ctest_ga3.c ga3_generated.h ga3.h
	gcc -I. -O3 ctest_ga3.c -o ctest_ga3 -lm

.PHONY:
clean:
	rm -f ctest_ga3 ga3_generated.h
