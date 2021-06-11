CC = gcc
CFLAGS = -O -std=gnu11 -fwrapv
LDLIBS = -lm

default: out.txt
	cat out.txt

out.txt: epsilon
	./epsilon > out.txt

epsilon: main.o tau.o
	$(CC) -o epsilon main.o tau.o $(LDLIBS)


main.o: main.c
	$(CC) $(CFLAGS) -c main.c 
tau.o: tau.c
	$(CC) $(CFLAGS) -c tau.c

.PHONEY: clean
clean: 
	$(RM) main.o tau.o epsilon out.txt

.PHONEY: test
test:
	echo $(LDLIBS)
	echo $(CC)
	echo $(RM)
