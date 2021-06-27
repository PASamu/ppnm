CC = gcc
CFLAGS = -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: out.txt
	cat out.txt

out.txt: main
	@./$< > $@

main: main.o roots.o qr.o
	$(CC) -o main main.o roots.o qr.o $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

roots.o: roots.c qr.h
	$(CC) $(CFLAGS) -c roots.c

.PHONEY: clean
	clean:
	$(RM) out.txt main *.o

