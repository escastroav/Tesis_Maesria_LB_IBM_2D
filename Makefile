CC = g++
CFLAGS = -I /usr/local/include/eigen3
DEPS = Waves_D2Q5.h IBM.h ComputeEpsilon.h
OBJ = main.o Waves_D2Q5.o IBM.o ComputeEpsilon.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main.x: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
