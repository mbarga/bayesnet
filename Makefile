CC = clang
CFLAGS = $(shell pkg-config --cflags glib-2.0 gtk+-2.0) -g -I./include -Wall -O0 
LDFLAGS = $(shell pkg-config --libs glib-2.0 gtk+-2.0) -lm
#CFLAGS = -g -I./include -Wall -O3 -std=c11 -pedantic  

all: main

main:
	$(CC) $(CFLAGS) ./src/readfile.c ./src/library.c ./src/probability.c ./src/score.c ./src/search.c ./src/main.c -o ./bin/hc $(LDFLAGS)
score:
	$(CC) $(CFLAGS) ./src/score.c ./src/test/bdetest.c -o ./bin/bdetest $(LDFLAGS)
prob:
	$(CC) $(CFLAGS) ./src/probability.c ./src/test/probtest.c -o ./bin/probtest $(LDFLAGS)
read:
	$(CC) $(CFLAGS) ./src/readfile.c ./src/test/testread.c -o ./bin/readtest $(LDFLAGS)

clean:
	rm -f *~.o ./bin/*
