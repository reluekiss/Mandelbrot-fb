CC = gcc
CFLAGS = -ggdb -static -O3 -Wall -Wextra
LDFLAGS = -lm -lpthread -ffunction-sections -fdata-sections -flto

all: clean build build/mandelbrot build/mandelbrot-arbitrary build/hopf-fibration

build/mandelbrot: src/mandelbrot.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

build/mandelbrot-arbitrary: src/mandelbrot-arbitrary.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) -lmpfr -lgmp

build/hopf-fibration: src/hopf-fibration.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

build:
	mkdir -p build 

clean:
	rm -rf build
