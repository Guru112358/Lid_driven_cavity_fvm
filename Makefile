CXX := g++
ACCFLAGS_1 := -O3 -march=native -fopt-info-vec


run: main.cpp
	${CXX} ${ACCFLAGS_1} -o run main.cpp 


clean:
	rm -f *.o *.csv *.png run *.plt *.dat
