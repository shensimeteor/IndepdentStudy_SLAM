
ALL: test.x test_obs.x Ceres_Burgers.x

Ceres_Burgers.x: Burgers.o Ceres_Burgers.o Observation.o
	g++ $^ -o $@

test.x: test.o Burgers.o
	g++ $^ -o $@

test_obs.x: test_obs.o Observation.o
	g++ $^ -o $@

test.o: Burgers.o
test_obs.o: Observation.o
Ceres_Burgers.o: Observation.o Burgers.o

%.o: %.cpp
	g++ -std=c++11 -c $<

.PHONY: clean
clean:
	rm -rf *.x *.o
