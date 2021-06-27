CFLAGS = -O -std=gnu11
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: fit.png

fit.png: fit.txt Makefile
	echo '\
	set terminal png;\
	set output $@;\
	set xlabel "x";\
	set ylabel "y";\
	set key outside;\
	set title "Neural network fit";\
	f(x)=x*exp(-x);\
	plot\
	 "$<" index 0 with points title "data"\
	,"$<" index 1 with lines title "fit"\
	'|pyxplot 

fit.txt: main
	@./$< > $@

main: main.o neurons.o qnewt.o



.PHONEY: clean
	clean: 
		$(RM) out.txt main *.png *.o

