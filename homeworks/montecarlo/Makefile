CC = gcc
CFLAGS = -Wall -O -std=gnu11
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: main
	@./$< > $@

main: main.o monte.o
	$(CC) -o $@ $^ $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -o $@ -c $<

monte.o: monte.c 
	$(CC) $(CFLAGS) -o $@ -c $<

.PHONEY: clean
	clean:
	$(RM) out.txt main *.o
