CC = g++
CFLAGS = 
DEPS = Waves_D2Q5.h IBM.h Rotor.h
OBJ = main.o Waves_D2Q5.o IBM.o Rotor.h

%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

main.x: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)
