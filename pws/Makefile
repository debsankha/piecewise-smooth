CCFLAGS=-I./ -O3 -I../numeric -ggdb

hardcol: ../numeric/librk4.o ../numeric/rk4.h hardcol.cpp libhardcol.o
	g++ $(CCFLAGS) -o hardcol.out hardcol.cpp ../numeric/librk4.o libhardcol.o 

libhardcol.o: libhardcol.cpp hardcol.h ../numeric/rk4.h
	g++ $(CCFLAGS) -c libhardcol.cpp 

probe-transient: probe_transient.cpp libhardcol.o ../numeric/librk4.o ../numeric/rk4.h
	g++ $(CCFLAGS) -o probe-transient.out probe_transient.cpp libhardcol.o ../numeric/librk4.o 

poincare: poincare.cpp 
	g++ $(CCFLAGS) -o poincare.out poincare.cpp ../numeric/librk4.o

shm-rk4: ../numeric/librk4.o shm-pw.o ../numeric/rk4.h
	g++ $(CCFLAGS) -o shm-rk4.out  ../numeric/librk4.o shm-pw.o

showdiff: ../numeric/librk4.o showdiff.cpp ../numeric/rk4.h
	g++ $(CCFLAGS) -o showdiff.out ../numeric/librk4.o showdiff.cpp 

testrk4: ../numeric/librk4.o testrk4.cpp ../numeric/rk4.h
	g++ $(CCFLAGS) -o testrk4.out ../numeric/librk4.o testrk4.cpp 

../numeric/librk4.o: ../numeric/librk4.cpp ../numeric/rk4.h
	g++ $(CCFLAGS) -c -o ../numeric/librk4.o  ../numeric/librk4.cpp

shm-pw.o: shm-pw.cpp
	g++ $(CCFLAGS) -c  shm-pw.cpp

scanvxarea:  scanvxarea.cpp
	g++ $(CCFLAGS) -o scanvxarea.out scanvxarea.cpp

clean:
	rm ../numeric/librk4.o shm-pw.o
