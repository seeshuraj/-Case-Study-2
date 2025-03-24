CXX = g++
CXXFLAGS = -std=c++11 -O3 -Wall -fopenmp
TARGET = gmres_solver

all: $(TARGET)

$(TARGET): gmres.cpp
	$(CXX) $(CXXFLAGS) -o $(TARGET) gmres.cpp

clean:
	rm -f $(TARGET) *.txt *.pdf

.PHONY: all clean