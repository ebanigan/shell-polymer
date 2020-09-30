all: runme

CPP = g++
OPTFLAGS = -O3
INCFLAGS= -I ~/boost_1_71_0

runme: main.cpp
	$(CPP) $(INCFLAGS) $(OPTFLAGS) -o $@ main.cpp
