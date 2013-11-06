TARGETS = goodsets goodsets128 goodsets128h


CFLAGS=-O3 -fomit-frame-pointer -funroll-all-loops -pthread -lrt
#CFLAGS=-g 

CC ?= gcc
#CC = /opt/open64/bin/opencc

all: $(TARGETS) 

goodsets: goodsets.c setops.o 
				 $(CC) $(CFLAGS) -o $@ $< setops.o -lpthread
#$(TARGETS): %: %.c setops.o 
#				 $(CC) $(CFLAGS) -o $@ $< setops.o -lpthread

goodsets128: goodsets128.c setops128.o 
				 $(CC) -DH128 $(CFLAGS) -o $@ goodsets128.c setops128.o -lpthread

goodsets128h: goodsets128.c setops.o 
				 $(CC) $(CFLAGS) -o $@ goodsets128.c setops.o -lpthread

setops128.o: setops128.c setops128.h
	 $(CC) $(CFLAGS) -c $<

setops.o: setops.c setops.h
	 $(CC) $(CFLAGS) -c $<

clean:
	rm -f $(TARGETS) *.o
