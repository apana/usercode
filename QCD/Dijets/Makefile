src     = ReadDijetTree
EXEC   = bin/$(src)
object = src/$(src).o  
source = src/$(src).cc
header = src/$(src).h

CXX           = g++
CXXFLAGS      = -O -Wall -fPIC
LD            = g++
ROOTCFLAGS   := $(shell root-config --cflags)
ROOTLIBS     := $(shell root-config --libs)
ROOTGLIBS    := $(shell $(ROOTSYS)/bin/root-config --glibs)
CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS) -lEG 
SHARE         = -shared


TARGET: $(EXEC)

$(EXEC): $(object)
	$(LD) -g -o $@ $(object) $(LIBS)

$(object): $(source) $(header)
	$(LD) -g -o $@ $(ROOTCFLAGS) -c $(source)

clean:
	rm $(object) $(EXEC)
