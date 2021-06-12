CC = gcc
CFLAGS = -O -Wall -std=gnu11
LDLIBS = -lm

.PHONEY:default
default: out.txt 
	@cat $<

out.txt: main
	@./main 1 2 3 4 5 > $@

main: main.o
	@$(CC) $^ -o $@ $(LDLIBS)

main.o: main.c
	@$(CC) $(CFLAGS) -c $<


.PHONEY: clean
clean:
	@$(RM) main 


