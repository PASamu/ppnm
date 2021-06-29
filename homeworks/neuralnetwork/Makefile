CC = gcc
CFLAGS = -O -std=gnu11 -Wall
CFLAGS += $(shell gsl-config --cflags) -pedantic
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: fit.png

fit.png: fit.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set xlabel "x";\
	set ylabel "y";\
	set key outside;\
	set title "Neural network fit";\
	plot\
	 "$<" index 0 using 1:2 with points title "data"\
	,"$<" index 1 using 1:2 with lines title "fit"\
	' | pyxplot 

fit.txt: main
	./$< > $@

main: main.o neurons.o qnewt.o neurons.h
	$(CC) $(CFLAGS) $^ -o $@ $(LDLIBS)

.PHONEY: clean
clean: 
	$(RM) main *.png *.o *.txt

