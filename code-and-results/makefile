
all: compile execute

compile:
	c++  -o main.out $(wildcard *.cpp) -larmadillo -std=c++11 -Ofast

execute:
	./main.out

plot1:
	python3 plot.py

perihelion:
	python3 perihelion.py

plotLE:
	python3 plotLE.py

plot3d:
	python3 plot3d.py

animation:
	python3 animate.py

get_angle:
	python3 get_angle.py

energy:
	python3 plotenergy.py
