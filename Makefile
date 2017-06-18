CC = g++
CFLAGS = -Ofast -funroll-loops -c
LFLAGS = -Ofast -funroll-loops
OMPFLAGS = -fopenmp
OBJS = LinearOpticalTransform.o main.o MeritFunction.o BFGS_Optimization.o

all: LinearOpticalSimulation

LinearOpticalSimulation: $(OBJS)
	$(CC) $(LFLAGS) $(OMPFLAGS) $(OBJS) -o LinearOpticalSimulation

main.o: main.cpp
	$(CC) $(CFLAGS) main.cpp

LinearOpticalTransform.o: LinearOpticalTransform.cpp
	$(CC) $(CFLAGS) $(OMPFLAGS) LinearOpticalTransform.cpp

MeritFunction.o: MeritFunction.cpp
	$(CC) $(CFLAGS) MeritFunction.cpp

BFGS_Optimization.o: BFGS_Optimization.cpp
	$(CC) $(CFLAGS) BFGS_Optimization.cpp

clean:
	rm *.o LinearOpticalSimulation
