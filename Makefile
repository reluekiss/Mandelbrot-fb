CC = gcc
CFLAGS = -ggdb -static -O3 -Wall -Wextra
LDFLAGS = -lm -lpthread -ffunction-sections -fdata-sections -flto

all: clean build build/mandelbrot build/mandelbrot-arbitrary build/hopf-fibration build/bott-periodicity build/framed-cobordism build/j-homomorphism build/torus-homology build/dec-mesh

build/mandelbrot: src/mandelbrot.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

build/mandelbrot-arbitrary: src/mandelbrot-arbitrary.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS) -lmpfr -lgmp

build/hopf-fibration: src/hopf-fibration.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

build/bott-periodicity: src/bott-periodicity.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

build/framed-cobordism: src/framed-cobordism.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

build/j-homomorphism: src/j-homomorphism.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

build/torus-homology: src/torus-homology.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

build/dec-mesh: src/dec-mesh.c
	$(CC) $(CFLAGS) $< -o $@ $(LDFLAGS)

build:
	mkdir -p build 

clean:
	rm -rf build
