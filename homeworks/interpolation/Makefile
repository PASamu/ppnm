CFLAGS = -O -Wall -Werror 
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)

default: linterp.png

linterp.png: out.txt
	echo '\
	set terminal png;\
	set output "$@";\
	set title "Linear spline";\
	set key top left;\
	set tics out;\
	set xlabel "x";\
	set ylabel "y";\
	plot \
	"points.txt" using 1:2 with points title "points",\
	,"$<" using 1:2 with line color green title "lin spline",\
	,"$<" using 1:3 with line color green title "integ lin spline",\
	,"$<" using 1:4 with line color red title "GSL lin spline",\
	,"$<" using 1:5 with line color red title "GSL integ lin spline",\
	' | pyxplot

out.txt: main
	./$< > $@

main: main.o linterp.o

.PHONEY: clean
clean:
	$(RM) *.txt *.o *.png
