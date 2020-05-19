# macros
CPP= g++
CFLAGS= -O3 -w -std=c++11 -fopenmp
EXECUTABLE=dtsppl
GRBPATH=/opt/gurobi901/linux64

# targets
all: data.o ilp_formulation_1.o ilp_formulation_2.o brkga.o main.o
	$(CPP) $(CFLAGS) -w -I$(GRBPATH)/include -L$(GRBPATH)/lib -o $(EXECUTABLE) data.o ilp_formulation_1.o ilp_formulation_2.o brkga.o main.o -lgurobi_g++5.2 -lgurobi90

data.o: data.cpp data.h
	$(CPP) $(CFLAGS) -c data.cpp

ilp_formulation_1.o: ilp_formulation_1.cpp ilp_formulation_1.h data.h
	$(CPP) $(CFLAGS) -I$(GRBPATH)/include -L$(GRBPATH)/lib -c ilp_formulation_1.cpp

ilp_formulation_2.o: ilp_formulation_2.cpp ilp_formulation_2.h data.h
	$(CPP) $(CFLAGS) -I$(GRBPATH)/include -L$(GRBPATH)/lib -c ilp_formulation_2.cpp

brkga.o: brkga.cpp brkga.h data.h
	$(CPP) $(CFLAGS) -I$(GRBPATH)/include -L$(GRBPATH)/lib -c brkga.cpp

main.o: main.cpp data.h ilp_formulation_1.h ilp_formulation_2.h brkga.h
	$(CPP) $(CFLAGS) -I$(GRBPATH)/include -L$(GRBPATH)/lib -c main.cpp

# remove
clean:
	rm -f ${EXECUTABLE} *.o
