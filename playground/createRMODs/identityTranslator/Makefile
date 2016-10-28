include ../../rose-hopper.mk
#include ../rose-carver.mk

#DBG         = -g
OPT         = -O3 
CFLAGS      = $(DBG) $(OPT)
CPP         = $(CXX)


EXTRA_LIBS  = 
LIBDIR      = .

SRCS        = \
	main.C \

OBJ_DIR	     =obj
SUBDIR	    = $(OBJ_DIR)/flops $(OBJ_DIR)/memory $(OBJ_DIR)/functions $(OBJ_DIR)/utils \
	$(OBJ_DIR)/loops


_OBJS	     =$(SRCS:%.C=%.o)
OBJS	     =$(patsubst %,$(OBJ_DIR)/%,$(_OBJS))

TARGET_BIN  = identityTrans

all: createOBJDIR $(TARGET_BIN) install

createOBJDIR	:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(SUBDIR)

$(OBJ_DIR)/%.o: %.C
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CPP) $(CFLAGS) $(CXXFLAGS)  -c -o $@ $< -I$(INCLUDE)
	@echo 'Finished building: $<'
	@echo ' '

$(TARGET_BIN) : $(OBJS)
	@echo 'Invoking C++ Linker'
	$(LIB_TOOL) --mode=link $(CPP) $(CFLAGS) $(CXXFLAGS) -o $@ $^ $(ROSE_LIBS)  -L$(LIBDIR) $(EXTRA_LIBS) 
	@echo 'Finished building target: $@'
	@echo ' '

install: $(TARGET_BIN)
	@echo $(TARGET_BIN) is installed
	@echo ' '
clean: 
	rm -f core *.o $(OBJS) *~  $(TARGET_BIN) 

clean-tests: 
	rm -f core *.o  *~ *.c.pdf *.c.dot rose_*.c rose_*.cpp mint_*.cpp mint_*.cu mint_*.c 

depend:
	$(CPP) -E -MM -I$(INCLUDE) $(EXTRA_CFLGS) $(SRCS) > DEPENDS

include DEPENDS
