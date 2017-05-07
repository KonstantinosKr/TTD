# C files
C_SRC=

# C++ files
CPP_SRC=main.cpp

# Program name
EXE=main

# Floating point type
REAL=double

# Debug version
DEBUG=no

TASK_CXX=tasksys.cpp
TASK_LIB=-lpthread
TASK_OBJ=tasksys.o
#TASK_OBJ=

ifeq ($(REAL),double)
  REAL=-DREAL=double
else
  REAL=-DREAL=float
endif

ifeq ($(DEBUG),yes)
  CFLAGS=-Iobjs/ -g -m64 $(REAL) -DDEBUG
else
  #CFLAGS=-Iobjs/ -O0 $(REAL) -no-vec
  CFLAGS=-Iobjs/ -O3 $(REAL) -xHOST -restrict -vec-report6 -std=c++0x  
endif

CC=icc $(CFLAGS)
CXX=icpc $(CFLAGS)

CPP_OBJS=$(addprefix objs/, $(CPP_SRC:.cpp=.o))
C_OBJS=$(addprefix objs/, $(C_SRC:.c=.o))
LIBS= -lpthread

default: dirs $(CPP_OBJS) $(C_OBJS) $(EXE)

.PHONY: dirs clean print

print:
	@echo $(CPP_OBJS)
	@echo $(C_OBJS)

dirs:
	/bin/mkdir -p objs/
	/bin/mkdir -p objs/krylov

del:
	find ./ -iname "*.dump" -exec rm '{}' ';'

clean:
	/bin/rm -rf objs *~ $(EXE)
	find ./ -iname "*.dump" -exec rm '{}' ';'

$(EXE): $(CPP_OBJS) $(C_OBJS) $(ISPC_OBJS)
	$(CXX) -o $@ $^ $(LIBS)

objs/%.o: %.cpp $(ISPC_HEADERS)
	$(CXX) $< -c -o $@

objs/%.o: %.c $(ISPC_HEADERS)
	$(CC) $< -c -o $@


