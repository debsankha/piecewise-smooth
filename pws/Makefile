shm-rk4: ../numeric/librk4.o shm-pw.o ../numeric/rk4.h
	g++ -ggdb -o shm-rk4.out  ../numeric/librk4.o shm-pw.o

../numeric/librk4.o: ../numeric/librk4.cpp
	g++ -ggdb -c -o ../numeric/librk4.o -I../numeric ../numeric/librk4.cpp

shm-pw.o: shm-pw.cpp
	g++ -ggdb -c -I../numeric shm-pw.cpp

clean:
	rm ../numeric/librk4.o shm-pw.o
