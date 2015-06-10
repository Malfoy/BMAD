#CC=/usr/bin/g++
CC=g++
CFLAGS=  -Wall  -O3 -std=c++11 -march=native -pthread
LDFLAGS=-pthread


ifeq ($(gprof),1)
CFLAGS=-std=c++0x -pg -O3  -march=native
LDFLAGS=-pg
endif

ifeq ($(valgrind),1)
CFLAGS=-std=c++0x -O3 -g
LDFLAGS=-g
endif


EXEC=minhashTest

all: $(EXEC)

minhashTest: minhashTest.o   minhash.o nw.o
	$(CC) -o $@ $^ $(LDFLAGS)

minhashTest.o: main.cpp minhash.h nw.h
	$(CC) -o $@ -c $< $(CFLAGS)

minhash.o: minhash.cpp minhash.h
	$(CC) -o $@ -c $< $(CFLAGS)

nw.o: nw.cpp nw.h
	$(CC) -o $@ -c $< $(CFLAGS)


clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
