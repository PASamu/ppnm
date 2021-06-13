CFLAGS = -Wall -pthread -Ofast -march=native
CC = gcc
LDLIBS = -lm

.PHONEY:default
default: out.txt
	cat out.txt

out.txt: main
	./$< > $@

.PHONEY:clean
clean: 
	$(RM) out* main

