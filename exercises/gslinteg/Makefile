CC = gcc
CFLAGS = -Wall -Werror -O1
CFLAGS += $(shell gsl-config --cflags)
LDLIBS += $(shell gsl-config --libs)

.PHONEY:default
default:out.txt
	cat out.txt

out.txt: main
	./$< > $@

main: main.o
	$(CC) -o main main.o $(LDLIBS)

main.o: main.c
	$(CC) $(CFLAGS) -c main.c

.PHONEY:clean
clean:
	$(RM) out* main 

.PHONEY: test
test:
	echo CFLAGS
	echo $(CFLAGS)
	echo $(LDLIBS)
