# C files
C_SRC=

# C++ files
CPP_SRC=gjk.cpp ttd.cpp

# Program name
EXE=ttd

# Floating point type
REAL=double

# Debug version
DEBUG=no

ifeq ($(REAL),double)
  REAL=-DiREAL=double
else
  REAL=-DiREAL=float
endif

ifeq ($(DEBUG),yes)
  CFLAGS=-Iobjs/ -g -m64 $(REAL) -DDEBUG
else
  #CFLAGS=-Iobjs/ -O0 $(REAL) -no-vec
  CFLAGS=-Iobjs/ -O3 $(REAL) -xHost -qopt-report5 -qopt-report-phase=vec
  #CFLAGS=-Iobjs/ -O3 $(REAL) -openmp -xHost -vec-report6 -override-limits 


  #CFLAGS=-Iobjs/ -O3 $(REAL) -ftree-vectorizer-verbose=2
  #CFLAGS=-Iobjs/ -O0 $(REAL) 
endif

CC=icc $(CFLAGS)
CXX=icpc $(CFLAGS)

CPP_OBJS=$(addprefix objs/, $(CPP_SRC:.cpp=.o))
C_OBJS=$(addprefix objs/, $(C_SRC:.c=.o))
LIBS= 

default: dirs $(ISPC_HEADERS) $(CPP_OBJS) $(C_OBJS) $(EXE)

.PHONY: dirs clean print

print:
	@echo $(CPP_OBJS)
	@echo $(C_OBJS)

dirs:
	/bin/mkdir -p objs/

del:
	find ./ -iname "*.dump" -exec rm '{}' ';'

clean:
	/bin/rm -rf objs *~ $(EXE)
	find ./ -iname "*.dump" -exec rm '{}' ';'

$(EXE): $(CPP_OBJS) $(C_OBJS)
	$(CXX) -o $@ $^ $(LIBS)

objs/%.o: %.cpp
	$(CXX) $< -c -o $@

