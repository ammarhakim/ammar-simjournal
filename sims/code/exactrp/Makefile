all: exacteulerrp reserv

exacteulerrp: exacteulerrp.o
	g++ exacteulerrp.o -o exacteulerrp

exacteulerrp.o: exacteulerrp.cc
	g++ -std=c++14 -O3 -Wall -c exacteulerrp.cc

reserv: reserv.o
	g++ reserv.o -o reserv

reserv.o: reserv.cc
	g++ -std=c++14 -O3 -Wall -c reserv.cc

clean:
	rm -rf exacteulerrp exacteulerrp.o reserv reserv.o
