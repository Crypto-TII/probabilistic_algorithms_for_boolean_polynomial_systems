
CC=gcc
CFLAGS= -Wall -pedantic
INC=    -I.
LIBS=   -lm
DEPS=   $(wildcard *.h) $(wildcard tests/*.h)
OBJ=    $(patsubst %.c,%.o,$(wildcard *.c)) $(patsubst %.c,%.o,$(wildcard tests/*.c))


%.o: %.c $(DEPS)
	$(CC) $(INC) -c -o $@ $< $(CFLAGS)

main: $(OBJ)
	$(CC) $(INC) -o $@ $^ $(CFLAGS) $(LIBS)

run: main
	./main

.PHONY: clean

clean:
	rm -f main *.o tests/*.o
