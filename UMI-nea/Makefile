CC = gcc
CXX=g++

UMI-nea: UMI-nea.main.cpp src/UMI-nea.o edlib/edlib.o
	$(CXX)  -std=c++14 -pthread -g -o UMI-nea UMI-nea.main.cpp src/UMI-nea.o edlib/edlib.o -I edlib/include/ -I include/

src/UMI-nea.o: src/UMI-nea.cpp include/UMI-nea.h
	cd src/ && $(CXX) -std=c++14  -c UMI-nea.cpp -I ../include -I ../edlib/include/

edlib/edlib.o:
	cd edlib/ && $(CXX) -c src/edlib.cpp -I include/

clean:
	rm -f UMI-nea edlib/edlib.o src/UMI-nea.o
