CC = clang
# -I for location of include files; -g for `gdb` debug; -ansi to adhere to c89 (ansi C) standards
CFLAGS = -g -I./include -Wall -O0 
#CFLAGS = -g -I./include -Wall -O3 -std=c11 -pedantic 

all: main

main:
# can list multiple src files :: $(CC) $(CFLAGS) ./src/template.c ./src/log.cc -o outname -lm
	$(CC) $(CFLAGS) ./src/library.c ./src/main.c -o ./bin/hc -lm

clean:
	rm -f *~.o ./bin/*

# for using object files and linking to personal shared library
#CFLAGS = -g -I./lib -L./lib -Wall -O3 -ansi -pedantic

#project: mylib.o source.o
#	$(CC) $(CFLAGS) -o project $^ -lm
#mylib.o: 
#	$(CC) $(CFLAGS) -c ./lib/mylib.c 
#source.o: 
#	$(CC) $(CFLAGS) -c ./src/source.c
#clean:
#	rm *.o project
