
all: compile execute

compile:
	c++ -Wall -o main.out $(wildcard *.cpp) -larmadillo -std=c++11

execute:
	./main.out

plot1:
	python3 plot.py
