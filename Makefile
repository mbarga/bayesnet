CC = clang
CFLAGS = $(shell pkg-config --cflags glib-2.0 gtk+-2.0) -g -I./include -Wall -O0 
LDFLAGS = $(shell pkg-config --libs glib-2.0 gtk+-2.0) -lm
#CFLAGS = -g -I./include -Wall -O3 -std=c11 -pedantic  

all: main

main:
	$(CC) $(CFLAGS) ./src/library.c ./src/main.c -o ./bin/hc $(LDFLAGS)

clean:
	rm -f *~.o ./bin/*