OBJECTS=main.o
CC=g++
INC_PATH=$(pwd)/include
CFLAGS=-Wall -Wextra -Werror -Wpedantic -Wno-unused-parameter -std=c++14 -O2 -I$(INC_PATH)
CLIBS= -L/usr/local/lib

all : main

main : $(OBJECTS)
	$(CC) $(CFLAGS) $^ $(CLIBS) -o $@

main.o : main.cpp
	$(CC) $(CFLAGS) -c $< -o $@


clean :
	rm -f *.o main

