CFLAGS = -O -std=gnu11
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
	f(x)=x*exp(-x);\
	plot\
	 "$(word 1, $<)" index 0 with points title "data"\
	,"$(word 1, $<)" index 1 with lines title "fit"\
	'|pyxplot 

fit.txt: main
	./$< > $@

main: main.o neurons.o qnewt.o


.PHONEY: clean
clean: 
	$(RM) main *.png *.o *.txt

