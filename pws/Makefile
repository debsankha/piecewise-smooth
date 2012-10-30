shm-rk4: ../numeric/librk4.o shm-pw.o ../numeric/rk4.h
	g++ -ggdb -o shm-rk4.out  ../numeric/librk4.o shm-pw.o

showdiff: ../numeric/librk4.o showdiff.cpp ../numeric/rk4.h
	g++ -ggdb -o showdiff.out ../numeric/librk4.o showdiff.cpp -I../numeric

testrk4: ../numeric/librk4.o testrk4.cpp ../numeric/rk4.h
	g++ -ggdb -o testrk4.out ../numeric/librk4.o testrk4.cpp -I../numeric

../numeric/librk4.o: ../numeric/librk4.cpp
	g++ -ggdb -c -o ../numeric/librk4.o -I../numeric ../numeric/librk4.cpp

shm-pw.o: shm-pw.cpp
	g++ -ggdb -c -I../numeric shm-pw.cpp

clean:
	rm ../numeric/librk4.o shm-pw.o
