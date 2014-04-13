DEBUG ?= 1
ifeq ($(DEBUG), 1)
	CFLAGS = -g3 -gdwarf2 -DDEBUG -DNPAR -I./include -Wall -O
	#CFLAGS = $(shell pkg-config --cflags glib-2.0 gtk+-2.0) -g3 -gdwarf2 -DDEBUG -I./include -Wall -O
	OMPCFLAGS = -fopenmp -std=c99 -g3 -DDEBUG -DPAR -I./include -Wall -O
else
	CFLAGS = -DNDEBUG -DNPAR -I./include -Wall -O3 -Wno-unknown-pragmas
	#CFLAGS = $(shell pkg-config --cflags glib-2.0 gtk+-2.0) -DNDEBUG -I./include -Wall -O3
	#CFLAGS = -g -I./include -Wall -O3 -std=c11 -pedantic
	OMPCFLAGS = -DNDEBUG -DPAR -fopenmp -std=c99 -g -I./include -Wall -O0
endif

LDFLAGS = $(shell pkg-config --libs glib-2.0 gtk+-2.0) -lm -lprofiler

CC = clang $(CFLAGS)
OMPCC = gcc $(OMPCFLAGS)

OMP_NUM_THREADS=4

all: main

main:
	$(CC) ./src/readfile.c ./src/BDE.c ./src/ran2.c ./src/util.c ./src/probability.c ./src/score.c ./src/search.c ./src/main.c -o ./bin/hc $(LDFLAGS)
omp:
	$(OMPCC) ./src/readfile.c ./src/BDE.c ./src/ran2.c ./src/util.c ./src/probability.c ./src/score.c ./src/search.c ./src/main.c -o ./bin/hc $(LDFLAGS)
	@echo "### OMP THREAD COUNT: ${OMP_NUM_THREADS} ###"
score:
	$(CC) ./src/score.c ./src/test/bdetest.c -o ./bin/bdetest $(LDFLAGS)
prob:
	$(CC) ./src/probability.c ./src/test/probtest.c -o ./bin/probtest $(LDFLAGS)
read:
	$(CC) ./src/readfile.c ./src/test/testread.c -o ./bin/readtest $(LDFLAGS)

clean:
	rm -f *~.o ./bin/*
