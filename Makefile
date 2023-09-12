CC = g++
CFLAGS = 
DEPS = Waves_D2Q5.h IBM.h
OBJ = main.o Waves_D2Q5.o IBM.o

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main.x: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
