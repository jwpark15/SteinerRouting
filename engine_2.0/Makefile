CC = g++
FLAGS= -std=c++11

engine: main.o MST.o PrimParams.o LRST.o KR.o results.o
	$(CC) $(FLAGS) main.o MST.o PrimParams.o LRST.o KR.o results.o -o engine

main.o: main.cpp MST.h LRST.h KR.h results.h
	$(CC) $(FLAGS) -c main.cpp

MST.o: MST.cpp MST.h PrimParams.h
	$(CC) $(FLAGS) -c MST.cpp

LRST.o: LRST.cpp LRST.h
	$(CC) $(FLAGS) -c LRST.cpp

KR.o: KR.cpp KR.h
	$(CC) $(FLAGS) -c KR.cpp

PrimParams.o: PrimParams.cpp PrimParams.h
	$(CC) $(FLAGS) -c PrimParams.cpp

results.o: results.cpp results.h
	$(CC) $(FLAGS) -c results.cpp

rand_benchmarks: BenchmarkGenerator.cpp
	$(CC) $(FLAGS) BenchmarkGenerator.cpp -o generate_benchmark

all: engine rand_benchmarks
