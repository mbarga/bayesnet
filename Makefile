CC = clang
CFLAGS = $(shell pkg-config --cflags glib-2.0 gtk+-2.0) -g -I./include -Wall -O0 
LDFLAGS = $(shell pkg-config --libs glib-2.0 gtk+-2.0) -lm
#CFLAGS = -g -I./include -Wall -O3 -std=c11 -pedantic  

all: main

main:
	$(CC) $(CFLAGS) ./src/readfile.c ./src/library.c ./src/bdescore.c ./src/main.c -o ./bin/hc $(LDFLAGS)
test:
	$(CC) $(CFLAGS) ./src/library.c ./src/bdescore.c ./src/main.c -o ./bin/hc $(LDFLAGS)
score:
	$(CC) $(CFLAGS) ./src/bdescore.c ./src/bdetest.c -o bdetest $(LDFLAGS)
prob:
	$(CC) $(CFLAGS) ./src/probability.c ./src/probtest.c -o probtest $(LDFLAGS)

clean:
	rm -f *~.o ./bin/*
