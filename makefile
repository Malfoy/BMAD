CC=/usr/bin/g++
#CC=g++
CFLAGS=  -Wall  -Ofast -std=c++11 -march=native -pthread
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

minhashTest: minhashTest.o  minhash.o nw.o utils.o xor.o
	$(CC) -o $@ $^ $(LDFLAGS)

minhashTest.o: main.cpp minhash.h nw.h utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

minhash.o: minhash.cpp minhash.h xor.h utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

utils.o: utils.cpp utils.h
	$(CC) -o $@ -c $< $(CFLAGS)

nw.o: nw.cpp nw.h
	$(CC) -o $@ -c $< $(CFLAGS)

xor.o: xor.cpp xor.h
	$(CC) -o $@ -c $< $(CFLAGS)

clean:
	rm -rf *.o
	rm -rf $(EXEC)


rebuild: clean $(EXEC)
