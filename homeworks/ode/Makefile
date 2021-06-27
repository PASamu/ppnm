CFLAGS = -Wall -O 
CFLAGS += $(shell gsl-config --cflags)
LDLIBS = -lm
LDLIBS += $(shell gsl-config --libs)


default: out.SIR.png 
	

out.SIR.png: sus.txt inf.txt rem.txt Makefile
	echo '\
	set terminal png;\
	set title "SIR-model";\
	set output "$@";\
	set key out;\
	set grid;\
	set xlabel "t (days)";\
	set ylabel "y (number of people)";\
	plot \
	 "sus.txt" using 1:2 with line title "S (Tc=2)"\
	,"sus.txt" using 1:3 with line title "S (Tc=4)"\
	,"sus.txt" using 1:4 with line title "S (Tc=7)"\
	,"inf.txt" using 1:2 with line title "I (Tc=2)"\
	,"inf.txt" using 1:3 with line title "I (Tc=4)"\
	,"inf.txt" using 1:4 with line title "I (Tc=7)"\
	,"rem.txt" using 1:2 with line title "R (Tc=2)"\
	,"rem.txt" using 1:3 with line title "R (Tc=4)"\
	,"rem.txt" using 1:4 with line title "R (Tc=7)"\
	' | pyxplot

sus.txt inf.txt rem.txt out.txt: main
	@./$< > $@

main: main.o ode.o
	$(CC) -o main main.o ode.o $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

ode.o: ode.c
	$(CC) $(CFLAGS) -c ode.c

.PHONEY: clean
clean:
	$(RM) main *.o *.txt *.png
