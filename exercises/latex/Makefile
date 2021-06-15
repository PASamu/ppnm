CFLAGS = -Wall -O1 -std=gnu11
LDLIBS = -lm

exp.pdf: exp.tex fig-pyxplot.pdf
	pdflatex $<
	pdflatex $<

out.data.txt:main
	./$< > $@

main:main.o exp.o

fig-pyxplot.pdf: fig-pyxplot.gpi out.data.txt
	pyxplot $<

clean:
	$(RM) exp.pdf *.log *.aux out* main *.o fig*.pdf fig*.tex
