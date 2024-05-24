all: lodepng.o main.o
	gcc *.o -o main

main.o: Основной_код.c lodepng.h
	gcc -c Основной_код.c

lodepng.o: lodepng.c lodepng.h
	gcc -c lodepng.c

clean:
	rm -f *.o
