CC = gcc
CFLAGS = -Wall -O -std=gnu11
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: main
	@./$< > $@

main: main.o integ.o
	$(CC) -o $@ $^ $(LDLIBS)

main.o: main.c integ.h
	$(CC) $(CFLAGS) -o $@ -c $<

integ.o: integ.c integ.h
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONEY: clean
	clean:
	$(RM) out.txt main *.o
