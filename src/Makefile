SRC=$(wildcard *.cpp)
HDR=$(wildcard *.h)
OBJ=bsvd.o binmat.o pbm.o GolombCoder.o coding.o util.o pnm.o #$(patsubst %.cpp,%.o,$(SRC))
TESTSRC=$(wildcard *_test.cpp *_tool.cpp)
TEST=$(patsubst %.cpp,%,$(TESTSRC))

CFLAGS=-O3 -Wall -fopenmp -msse4 -Wextra
#CFLAGS=-O -Wall -g -fopenmp

all: $(OBJ) $(TEST)

%.o: %.cpp $(HDR) Makefile
	echo $(HDR)
	g++ $(CFLAGS) -c -o $@ $<

%_test: %_test.cpp $(HDR) $(OBJ) Makefile
	g++  $(CFLAGS) -Wall -o $@ $< $(OBJ) -lgsl -lgslcblas -lm

%_tool: %_tool.cpp $(HDR) $(OBJ) Makefile
	g++  $(CFLAGS) -Wall -o $@ $< $(OBJ) -lgsl -lgslcblas -lm

clean:
	-rm $(TEST)
	-rm *.o
