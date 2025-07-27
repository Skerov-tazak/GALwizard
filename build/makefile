result: matrix.o functions.o testing.o 
	g++ -g matrix.o functions.o testing.o -o result 

matrix.o: matrix.cpp
	g++ -g -c matrix.cpp -o matrix.o

functions.o: functions.cpp
	g++ -g -c functions.cpp -o functions.o

testing.o: testing.cpp
	g++ -g -c testing.cpp -o testing.o

clean:
	rm matrix.o
	rm functions.o
	rm testing.o
	rm result
