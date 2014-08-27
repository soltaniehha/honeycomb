CXX = pgCC
CXX = g++32
CXX = g++
RM = rm -f
LDFLAGS = -DWITH_LAPACK -lblas -llapack
INCLUDES = -D__USE_GNU -I.. -I. -I/Users/mohammad/include
CXXCOMPILE = $(CXX) -DWITH_LAPACK $(INCLUDES)

SRCS = simpleHexagonal.cc
OBJS = $(subst .cc,.o,$(SRCS))

all: $(OBJS)
	$(CXX) $(LDFLAGS) -o  run $(OBJS) $(INCLUDES)

%.o: %.cc
	$(CXXCOMPILE) -g -c $<

clean:
	$(RM) $(OBJS)
