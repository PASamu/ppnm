CC = gcc
CFLAGS = -Wall -std=gnu11 -DNDEBUG
CFLAGS += $(shell /usr/bin/gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell /usr/bin/gsl-config --libs)

default: out.txt
	cat out.txt

out.txt: main
	./$< > $@

main: main.c

clean:
	$(RM) main out.txt
