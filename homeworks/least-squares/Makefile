CC = gcc
CFLAGS += -O -march=native -Wall -std=gnu11 $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

default: data.txt fitplot.pyxplot.png

fitplot.pyxplot.png: data.txt Makefile
	echo '\
	set terminal png;\
	set output "$@";\
	set key top right;\
	set tics out;\
	set xlabel("time (days)");\
	set ylabel("ln(activity)");\
	set title "Least square fit";\
	plot \
	  "$<" index 0 with errorbars title "Data points"\
	 ,"$<" index 1 using 1:2 with line title "Least square fit"\
	 ,"$<" index 1 using 1:3 with line title "Fitplus"\
	 ,"$<" index 1 using 1:4 with line title "Fitminus"\
	 '| pyxplot
	
data.txt: main
	./$< 1>$@ 2>log

main: main.o QR_decomp.o leastsq.o

QR_decomp.o: QR_decomp.c
	$(CC) $(CFLAGS) -c -o $@ $<

.PHONEY:
	clean: $(RM) main *.o plottingdata.txt out.txt measurements.txt
