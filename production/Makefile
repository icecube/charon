INCDIR = -I$(shell pythia8-config --includedir)
LIBS = $(shell pythia8-config --libs) $(shell gsl-config --libs)
CC = g++
CFLAGS = -std=c++11 -O2 -pedantic -W -Wall -Wshadow -fPIC -Wl,-rpath, $(shell gsl-config --cflags)
SRC = $(wildcard ./*.cc)
TARGET = $(patsubst %.cc,%,$(SRC))
OUT = $(patsubst %.cc,%.out,$(SRC))
# HEADERS = $(wildcard $(INCDIR)/*.h)
.PHONY: all clean run
all: $(TARGET)    
	@echo "DONE!"
$(TARGET): $(SRC)
	$(CC) $(CFLAGS) $@.cc -o $@ $(INCDIR) $(LIBS)
clean:
	rm -rf $(TARGET) 
