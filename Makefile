#desactivation of the implicite rules
.SUFFIXES:

#Definitions
EXECUTABLE = ECOGEN
CXX = mpic++
CXXFLAGS = -O3 -std=c++11
#CXXFLAGS = -g -std=c++11
# LDFLAGS =

dirs = $(shell find . -type d)
SOURCES = $(foreach dir,$(dirs),$(wildcard $(dir)/*.cpp))
OBJETS = $(SOURCES:.cpp=.o)

all: $(OBJETS)
		$(CXX) $^ -o $(EXECUTABLE) $(CXXFLAGS)

%o: %cpp
		$(CXX) -c $< -o $@ $(CXXFLAGS)


###

depend:
		makedepend $(SOURCES)

clean:
		rm -rf $(OBJETS)

cleanres:
		rm -rf ./results/*

#Creation of the executable


# DO NOT DELETE


