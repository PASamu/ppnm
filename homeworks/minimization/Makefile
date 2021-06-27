CC = gcc
CFLAGS = -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: out.txt
	cat out.txt

out.txt: main
	@./$< > $@

main: main.o minimization.o
	$(CC) -o main main.o minimization.o $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

minimization.o: minimization.c
	$(CC) $(CFLAGS) -c minimization.c

.PHONEY: clean
	clean:
	$(RM) out.txt main *.o
