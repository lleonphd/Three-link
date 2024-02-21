all: main.o
	gcc -o main main.o -lm

main.o: main.c
	gcc -c main.c

clean:
	rm main main.o
