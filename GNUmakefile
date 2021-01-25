CXX = g++
FCC = gfortran
LIBS = 

# Directory definitions
INCDIR = include
SRCDIR = src
OBJDIR = obj
EXECDIR = .
FSRCDIR = talysInterface/source
FOBJDIR = talysInterface/obj
TARGET = newtonRun

#Compiler flags
# Enable Root Usage, Include Headers, Usage of Fortran Code in C++
CXXFLAGS = `root-config --libs --cflags` -I$(INCDIR) -lgfortran -no-pie -fPIC -fpermissive
FCCFLAGS = -I$(TALYS)/source/

DEPS = xscn.h flux.h detector.h event.h exstates.h global.h organizer.h $(TALYS)/source/talys.cmb

SRC = newtonRun.cc xscn.cc flux.cc detector.cc event.cc exstates.cc organizer.cc

# CPP and FORTRAN OBJECTS
OBJS = $(SRC:%.cc=$(OBJDIR)/%.o)
FOBJS = $(shell find $(FOBJDIR) -name '*.o')

# TALYS objects, compiled when installing TALYS
TALYSO = $(filter-out $(TALYS)/source/talys.o, $(wildcard $(TALYS)/source/*.o))

all: $(TARGET)

# Compile c++ objects
$(OBJDIR)/%.o : $(SRCDIR)/%.cc
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

# Compile fortran objects (used to interface with TALYS)
$(FOBJDIR)/%.o : $(FSRCDIR)/%.f
	$(FCC) -c -o $@ $^ $(FCCFLAGS)

# Make target by linking it to fortran, c++ and TALYS objects 
$(TARGET) : $(OBJS) $(FOBJS)
	$(CXX) $(TALYSO) $^ -o $@ $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o $(FOBJDIR)/*.o *~ core 
