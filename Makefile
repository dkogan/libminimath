CC=gcc-4.5
CFLAGS = -Wall -Wextra -std=gnu99 -ffast-math -O3 -mtune=core2


minimath_generated.h: minimath_generate.pl
	./minimath.pl
