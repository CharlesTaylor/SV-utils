CXX ?= g++
OPT=-O2 -g
LDFLAGS += -lm 
CFLAGS +=  -Wall -Wextra -std=c++14
SOURCES= vcf2bed.cpp
OBJECTS=  $(SOURCES:.cpp=.o)
EXEC = vcf-to-bed


$(EXEC): $(OBJECTS)
	$(CXX) $(OBJECTS) $(OPT) -o $@ $(LDFLAGS)

.c.o:
	$(CXX) -c $(CFLAGS) $< -o $@

