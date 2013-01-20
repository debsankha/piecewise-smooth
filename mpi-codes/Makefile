CC=mpicc	#compile at any cluster node
CCFLAGS=-I./ -O3 -I../numeric

hardcol: ../numeric/librk4.o ../numeric/rk4.h hardcol.cpp libhardcol.o
	CC $(CCFLAGS) -o hardcol.out hardcol.cpp ../numeric/librk4.o libhardcol.o 

libhardcol.o: libhardcol.cpp hardcol.h ../numeric/rk4.h
	CC $(CCFLAGS) -c libhardcol.cpp 

probe-transient: probe_transient.cpp libhardcol.o ../numeric/librk4.o ../numeric/rk4.h
	CC $(CCFLAGS) -o probe-transient.out probe_transient.cpp libhardcol.o ../numeric/librk4.o 

poincare: poincare.cpp 
	CC $(CCFLAGS) -o poincare.out poincare.cpp ../numeric/librk4.o

shm-rk4: ../numeric/librk4.o shm-pw.o ../numeric/rk4.h
	CC $(CCFLAGS) -o shm-rk4.out  ../numeric/librk4.o shm-pw.o

showdiff: ../numeric/librk4.o showdiff.cpp ../numeric/rk4.h
	CC $(CCFLAGS) -o showdiff.out ../numeric/librk4.o showdiff.cpp 

testrk4: ../numeric/librk4.o testrk4.cpp ../numeric/rk4.h
	CC $(CCFLAGS) -o testrk4.out ../numeric/librk4.o testrk4.cpp 

../numeric/librk4.o: ../numeric/librk4.cpp
	CC $(CCFLAGS) -c -o ../numeric/librk4.o  ../numeric/librk4.cpp

shm-pw.o: shm-pw.cpp
	CC $(CCFLAGS) -c  shm-pw.cpp

clean:
	rm ../numeric/librk4.o shm-pw.o
