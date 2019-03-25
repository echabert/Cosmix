ObjSuf        = o
SrcSuf        = cc
ExeSuf        = run
LogSuf        = log

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs) 

CXX           = g++
#CXXFLAGS      = -std=c++0x
CXXFLAGS     += -O -Wall -fPIC $(DEFINES) -Wno-unused-result -Wshadow
CXXFLAGS     += $(ROOTCFLAGS) -I./

LD            = g++ 
LDFLAGS       = -g -O -Wall -fPIC -Wl,--no-undefined 

SOFLAGS       = -shared
LIBS          =  

#------------------------------------------------------------------------------

SOURCES       = $(wildcard *.$(SrcSuf))
OBJECTS       = $(SOURCES:.$(SrcSuf)=.$(ObjSuf))
EXECUTABLES   = $(SOURCES:.$(SrcSuf)=.$(ExeSuf))
LOGS          = $(SOURCES:.$(SrcSuf)=.$(LogSuf))

#------------------------------------------------------------------------------




libCosmix.so: DataAnalysis.cc CosmixSimu.cc DataAnalysis.h CosmixSimu.h
	g++ $(CXXFLAGS) $(ROOTLIBS) -shared -fPIC DataAnalysis.cc CosmixSimu.cc -o libCosmix.so 

SimuMain: libCosmix.so SimuMain.cc
	g++ $(CXXFLAGS) $(ROOTLIBS) -L`pwd` -lCosmix SimuMain.cc -o SimuMain

AnaMain: libCosmix.so AnaMain.cc
	g++ $(CXXFLAGS) $(ROOTLIBS) -L`pwd` -lCosmix AnaMain.cc -o AnaMain

clean:
	@echo "Cleaning..."
	@rm -f *.$(ObjSuf) *.$(ExeSuf) *.$(LogSuf)

#------------------------------------------------------------------------------


%.$(ExeSuf): %.$(SrcSuf) ../../.vectorDictionnary_C.so
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LIBS) $(GCCPARSER)


