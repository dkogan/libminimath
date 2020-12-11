TARGET = minimath_generated.h
HEADERS = $(TARGET) minimath.h

all: $(TARGET)

$(TARGET): minimath_generate.pl
	./$< > $@.tmp && mv $@.tmp $@

unittest: unittest.o
unittest.o: $(HEADER)
CFLAGS = -Wall -Wextra -std=gnu99 -ffast-math -O3 -I.


ifdef DESTDIR
install: $(TARGET)
	mkdir -p $(DESTDIR)/usr/include/
	install -m 0644 $(HEADERS) $(DESTDIR)/usr/include/
else
install:
	@echo "make install is here ONLY for the debian package. Do NOT run it yourself" && false
endif

check: unittest
	./$<

clean:
	rm -rf unittest unittest.o $(TARGET)

.PNONY: clean install check
