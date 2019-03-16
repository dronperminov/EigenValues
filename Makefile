COMPILER=g++
FLAGS=-Wall -O3
SRC=Matrix.cpp Vector.cpp main.cpp
TARGET=eigenvalues

all:
	$(COMPILER) $(FLAGS) $(SRC) -o $(TARGET)

clean:
	rm -rf $(TARGET)