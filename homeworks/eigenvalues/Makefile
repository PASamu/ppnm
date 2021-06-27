CC = gcc
CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: out.txt
	cat out.txt

out.txt: main
	./$< > $@

main: main.o jacobi.o
	$(CC) -o main main.o jacobi.o $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

jacobi.o: jacobi.c
	$(CC) $(CFLAGS) -c jacobi.c

clean:
	$(RM) *.o out.txt
