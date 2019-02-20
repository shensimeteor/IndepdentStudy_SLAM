
ALL: test.x test_obs.x run_burgers.x

FCFLAGS=-g -std=c++11 
LDFLAGS=-L/usr/local/lib -lconfig++

#Ceres_Burgers.x: Burgers.o Ceres_Burgers.o Observation.o
#	g++ $^ $(LDFLAGS) -o $@

test.x: test.o Burgers.o
	g++ $^ $(LDFLAGS) -o $@

test_obs.x: test_obs.o Observation.o
	g++ $^ $(LDFLAGS) -o $@

run_burgers.x: run_burgers.o Burgers.o
	g++ $^ $(LDFLAGS) -o $@

test.o: Burgers.o
test_obs.o: Observation.o
Ceres_Burgers.o: Observation.o Burgers.o
run_burgers.o: Burgers.o

%.o: ../src/%.cpp
	g++ $(FCFLAGS) -c $<

.PHONY: clean
clean:
	rm -rf *.x *.o
