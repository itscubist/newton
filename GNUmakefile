CXX = g++
LIBS = 

INCDIR = include
SRCDIR = src
OBJDIR = obj
EXECDIR = .
TARGET = newtonTest

# Enable Root Usage and Include
CXXFLAGS = `root-config --libs --cflags` -I$(INCDIR)

#DEPS = $(find $(INCDIR) -name '*.h')
DEPS = xscn.h flux.h detector.h event.h global.h

#SRC = $(find $(SRCDIR)/ -name '*.cc')
SRC = newtonTest.cc xscn.cc flux.cc detector.cc event.cc

OBJS = $(SRC:%.cc=$(OBJDIR)/%.o)

all: $(TARGET)

$(OBJDIR)/%.o : $(SRCDIR)/%.cc
	$(CXX) -c -o $@ $^ $(CXXFLAGS)

$(TARGET) : $(OBJS) 
	$(CXX) $^ -o $@ $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o *~ core 
