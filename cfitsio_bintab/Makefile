CC=gcc

CFLAGS=-I/arxiv/libraries/ubuntu/gcc/cfitsio3370/include/ -std=c99 -pedantic -Wall -Wextra -g
LDFALGS=-L//arxiv/libraries/ubuntu/gcc/cfitsio3370/lib/
LIBS=-lcfitsio -lm

default:write_bintab.exe

%.exe:%.o
	$(CC) $(LDFALGS) -o $@ $^ $(LIBS)
%.o:%.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm *.o *.exe
