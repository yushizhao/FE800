#
# created by Bartek 'BaSz' Szurgot (http://www.baszerr.eu)
#

CXX:=g++
CXXFLAGS:=-std=c++11

#CXX:=g++-4.6
#CXXFLAGS:=-std=c++0x

CXXFLAGS+=-O3 -DNDEBUG
#CXXFLAGS+=-g3


# common part
DOXYGEN:=doxygen
STRIP:=strip
CXXFLAGS+=-Wall -pipe -MD
LDFLAGS:=-lPocoNet -lPocoUtil -lPocoFoundation

SRCS:=$(wildcard *.cpp)
OBJS:=$(SRCS:.cpp=.o)
DEPS:=$(SRCS:.cpp=.d)

.PHONY: all
all: example.out

.PHONY: doc
doc: Doxyfile
	$(DOXYGEN) > /dev/null

Doxyfile:
	$(DOXYGEN) -g > /dev/null

-include $(DEPS)

%.o: %.cpp Makefile
	$(CXX) $(CXXFLAGS) -c $*.cpp

example.out: $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)
	$(STRIP) $@

.PHONY: clean
clean:
	rm -rfv html latex Doxyfile
	rm -vf *.o *.out *.d
