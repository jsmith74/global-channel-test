CC = g++
CFLAGS = -Ofast -funroll-loops -c
LFLAGS = -Ofast -funroll-loops
OMPFLAGS = -fopenmp
OBJS = LinearOpticalTransform.o main.o MeritFunction.o BFGS_Optimization.o

all: LinearOpticalSimulation Script DataProcess

DataProcess: dataProcess.cpp
	$(CC) dataProcess.cpp -o DataProcess

Script: script.cpp
	$(CC) $(OMPFLAGS) script.cpp -o Script

LinearOpticalSimulation: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o LinearOpticalSimulation

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

LinearOpticalTransform.o: LinearOpticalTransform.cpp
	$(CC) $(CFLAGS) LinearOpticalTransform.cpp

MeritFunction.o: MeritFunction.cpp
	$(CC) $(CFLAGS) MeritFunction.cpp

BFGS_Optimization.o: BFGS_Optimization.cpp
	$(CC) $(CFLAGS) BFGS_Optimization.cpp

clean:
	rm *.o LinearOpticalSimulation *.dat Script DataProcess
