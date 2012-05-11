CXX = g++
CFLAGS = -O3 -Wno-deprecated -ldai
LDFLAGS =
OBJS = FactoredMDP.o XML.o Experiment.o

OPTS = $(OBJS) $(CFLAGS) $(LDFLAGS)

default: depend main

$(OBJS): %.o : %.cpp
	$(CXX) $(CFLAGS) -o $@ -c $<

main: $(OBJS) main.cpp
	$(CXX) -o main main.cpp $(OPTS)

depend:
	$(CXX) $(CFLAGS) -MM *.cpp > Makefile.dep

clean:
	rm -f *.o main

include Makefile.dep
