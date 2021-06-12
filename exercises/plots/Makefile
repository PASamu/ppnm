#CFLAGS += `gsl-config --cflags`
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = $(shell gsl-config --libs)

.PHONEY: default
default: erf.pyxplot.png gamma.pyxplot.png

gamma.pyxplot.png: gamma.txt Makefile tab.gamma.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics in;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Gamma function";\
	plot \
	 "$<" using 1:2 with line title "tgamma from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-gamma"\
	,"$<" using 1:4 with line title "gamma approx"\
	,"tab.gamma.txt" using 1:2 with points title "tabulated data"\
	' | tee tmp2.gpi | pyxplot

erf.pyxplot.png: erf.txt Makefile tab.erf.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set key top left;\
	set tics in;\
	set xlabel "x";\
	set ylabel "y";\
	set title "Error function";\
	plot \
	 "$<" using 1:2 with line title "erf from math.h"\
	,"$<" using 1:3 with line title "gsl-sf-erf"\
	,"$<" using 1:4 with line title "erf approx"\
	,"tab.erf.txt" using 1:2 with points title "tabulated data"\
	' | tee tmp.gpi | pyxplot

tmp2.gpi: gamma.pyxplot.png

gamma.txt: gamma
	./$< > $@

gamma: gamma.o

gamma.o: gamma.c

tmp.gpi: erf.pyxplot.png

erf.txt: erf
	./$< > $@

erf: erf.o

erf.o: erf.c

.PHONEY:clean
clean: 
	$(RM) *.o erf gamma *.png tmp.gpi tmp2.gpi

.PHONEY:test
test:
	echo $(CFLAGS)
	echo $(LDLIBS)
