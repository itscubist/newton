CXX = g++
LIBS = 

INCDIR = include
SRCDIR = src
OBJDIR = obj
EXECDIR = .
TARGET = newtonTest

CXXFLAGS = `root-config --libs --cflags` -I$(INCDIR)

DEPS = $(find $(INCDIR) -name '*.h')
 
#SRC = $(find $(SRCDIR)/ -name '*.cc')
SRC = newtonTest.cc

OBJS = $(SRC:%.cc=$(OBJDIR)/%.o)

all: $(TARGET)

$(OBJDIR)/%.o : $(SRCDIR)/%.cc $(DEPS)
	$(CXX) -c -o $@ $< $(CXXFLAGS)

$(TARGET) : $(OBJS) 
	$(CXX) $^ -o $@ $(CXXFLAGS)

.PHONY: clean

clean:
	rm -f $(OBJDIR)/*.o *~ core $(INCDIR)/*~ 
