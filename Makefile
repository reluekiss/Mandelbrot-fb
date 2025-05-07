CC = gcc
CFLAGS = -static -O3 -Wall -Wextra
LDFLAGS = -lm -lpthread -ffunction-sections -fdata-sections -flto

all: clean build build/mandelbrot build/mandelbrot-arbitrary run

run:
	./build/mandelbrot

build/mandelbrot: src/main.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

build/mandelbrot-arbitrary: src/arbitrary.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) -lmpfr -lgmp

build:
	mkdir -p build 

clean:
	rm -rf build
