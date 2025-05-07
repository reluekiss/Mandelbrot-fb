CC = gcc
CFLAGS = -static -O3 -Wall -Wextra
LDFLAGS = -lm -lpthread -ffunction-sections -fdata-sections -flto

all: clean build build/mandelbrot run

run:
	./build/mandelbrot

build/mandelbrot: src/main.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

build:
	mkdir -p build 

clean:
	rm -rf build
