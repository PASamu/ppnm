CFLAGS = -Wall -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)
CC = gcc

default: out.diffeq.png

out.diffeq.png: out.txt Makefile
	echo '\
	set terminal png;\
	set title "Plot of Re(y) and Im(y)";\
	set output "$@";\
	set key out;\
	set grid;\
	set xlabel "t";\
	set ylabel "y";\
	plot \
	 "out.txt" using 1:3 with line title "Re(y(Re(t)))"\
	,"out.txt" using 2:3 with line title "Re(y(Im(t)))"\
	,"out.txt" using 2:4 with line title "Im(y(Im(t)))"\
	,"out.txt" using 1:4 with line title "Im(y(Re(t)))"\
	' | pyxplot

	
out.txt: main
	./$< > $@

main: main.o ode.o

.PHONEY: clean
clean: 
	$(RM) main *.o *.txt
