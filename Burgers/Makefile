
test.x: test.o Burgers.o
	g++ $^ -o $@

test.o: test.cpp
Burgers.o: Burgers.cpp

%.o: %.cpp
	g++ -c $<
