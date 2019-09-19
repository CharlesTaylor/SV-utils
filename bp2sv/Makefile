CXX ?= g++
OPT=-O2 -g
LDFLAGS += -lm -std=c++14
CFLAGS +=  -Wall -Wextra 
SOURCES= bp2sv.cpp
OBJECTS=  $(SOURCES:.cpp=.o)
EXEC = bp-to-sv


$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) $(OPT)  -o $@ $(LDFLAGS)

.c.o:
	$(CXX) -c $(CFLAGS)   $< -o $@

