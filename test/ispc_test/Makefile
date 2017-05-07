# C files
C_SRC=

# C++ files
CPP_SRC=main.cpp

# ISPC files
ISPC_SRC= test.ispc

# ISPC targets
ISPC_TARGETS=avx

# Program name
EXE=ttd

# Floating point type
REAL=double

# Debug version
DEBUG=no


ifeq ($(REAL),double)
  REAL=-DREAL=double
else
  REAL=-DREAL=float
endif

ifeq ($(DEBUG),yes)
  CFLAGS=-Iobjs/ -g -m64 $(REAL) -DDEBUG
  ISPC=ispc -g --arch=x86-64 $(REAL) -DDEBUG
else
  CFLAGS=-Iobjs/ -O3 -m64 $(REAL)
  ISPC=ispc -O3 --arch=x86-64 --woff $(REAL) 
endif

CC=icc $(CFLAGS)
CXX=icpc $(CFLAGS)

ISPC_OBJS=$(addprefix objs/, $(ISPC_SRC:.ispc=_ispc.o))
ISPC_HEADERS=$(addprefix objs/, $(ISPC_SRC:.ispc=_ispc.h))
CPP_OBJS=$(addprefix objs/, $(CPP_SRC:.cpp=.o))
C_OBJS=$(addprefix objs/, $(C_SRC:.c=.o))
LIBS= -lpthread

default: dirs $(ISPC_HEADERS) $(CPP_OBJS) $(C_OBJS) $(EXE)

.PHONY: dirs clean print

print:
	@echo $(ISPC_HEADERS)
	@echo $(CPP_OBJS)
	@echo $(C_OBJS)
	@echo $(ISPC_OBJS)

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

objs/%_ispc.h objs/%_ispc.o: %.ispc
	$(ISPC) --target=$(ISPC_TARGETS) $< -o objs/$*_ispc.o -h objs/$*_ispc.h
