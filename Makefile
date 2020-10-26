CC = g++ -std=c++17 -Wall # -g
#CC = g++ -pg
#CC = g++ -Wall 
CCFLAGS=  #-m64

CFLAGS=-Iparser/include -Lparser

PROG = metacount
SOURCES= utilities.cpp options.cpp  metacount.cpp core.cpp parser.cpp types.cpp
OBJECTS= $(SOURCES:.cpp=.o)
HEADERS= $(SOURCES:.cpp=.h)

all: $(PROG)

%.o:%.cpp   $(SOURCES) types.h
	$(CC) $(CCFLAGS) $(CFLAGS) $< -c -o $@  

clean:
	rm -rf $(OBJECTS) $(PROG)

$(PROG): $(OBJECTS) $(HEADERS) types.h
	$(CC) $(CCFLAGS) $(CFLAGS)  $(OBJECTS) -o $(PROG) -lbamparser -lz -lstdc++fs

