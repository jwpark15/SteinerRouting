CC = g++
FLAGS= -std=c++11

engine: main.o MST.o PrimParams.o LRST.o KR.o results.o
	$(CC) main.o MST.o PrimParams.o LRST.o KR.o results.o -o engine

main.o: main.cpp MST.h LRST.h KR.h results.h
	$(CC) -c main.cpp

MST.o: MST.cpp MST.h PrimParams.h
	$(CC) -c MST.cpp

LRST.o: LRST.cpp LRST.h
	$(CC) -c LRST.cpp

KR.o: KR.cpp KR.h
	$(CC) -c KR.cpp

PrimParams.o: PrimParams.cpp PrimParams.h
	$(CC) -c PrimParams.cpp

results.o: results.cpp results.h
	$(CC) -c results.cpp
