CC=gcc-4.5
CFLAGS = -MMD -Wall -Wextra -std=gnu99 -ffast-math -O3 -mtune=core2

all: minimath_generated.h

minimath_generated.h: minimath_generate.pl
	./minimath_generate.pl

unittest: unittest.o

clean:
	rm -rf *.o unittest

.PNONY: clean

-include *.d
