include ../rose.mk

#DBG         = -g
OPT         = -O3 
CFLAGS      = $(DBG) $(OPT)
CPP         = $(CXX)


EXTRA_LIBS  = 
LIBDIR      = . -L/opt/boost/compileTree/boost-1-50-0/lib/ -lboost_iostreams -lboost_system -lboost_python -L/usr/lib64/python2.6/config/ -lpython2.6

SRCS        = \
	main.C \
	./scalars/Scalars.C \
	./flops/Flops.C \
	./memory/ReadWriteProp.C \
	./memory/AccessAnalysis.C \
	./memory/FuncSideEffect.C \
	./functions/FunctionHandler.C \
	./handler/ScopeStmtHandler.C \
	./functions/FunctionCallHandler.C \
	./loops/LoopHandler.C \
	./conditionals/IfHandler.C \
	./specialBasicBlocks/specialBBHandler.C\
	./modules/ModuleHandler.C \
	./wholearrays/WholeArrayExp.C \
	./utils/Utils.C \
	./utils/ArrayUtils.C \
	./utils/OutputFormatting.C \
	./utils/ExaSatOptions.C \
	./DAG/DAG.C \

OBJ_DIR	     =obj
SUBDIR	    = $(OBJ_DIR)/scalars $(OBJ_DIR)/flops $(OBJ_DIR)/memory $(OBJ_DIR)/handler $(OBJ_DIR)/functions $(OBJ_DIR)/utils \
	$(OBJ_DIR)/loops $(OBJ_DIR)/specialBasicBlocks $(OBJ_DIR)/modules $(OBJ_DIR)/wholearrays $(OBJ_DIR)/conditionals $(OBJ_DIR)/DAG



_OBJS	     =$(SRCS:%.C=%.o)
OBJS	     =$(patsubst %,$(OBJ_DIR)/%,$(_OBJS))

TARGET_BIN  = genCodeDescript

all: createOBJDIR $(TARGET_BIN) install

createOBJDIR	:
	mkdir -p $(OBJ_DIR)
	mkdir -p $(SUBDIR)

$(OBJ_DIR)/%.o: %.C
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CPP) $(CFLAGS) $(CXXFLAGS)  -c -o $@ $< $(ROSE_INCLUDE)
	@echo 'Finished building: $<'
	@echo ' '

$(TARGET_BIN) : $(OBJS)
	@echo 'Invoking C++ Linker'
	$(LIB_TOOL) --mode=link --tag=LD $(CPP) $(CFLAGS) $(CXXFLAGS) -o $@ $^ $(ROSE_LIBS) $(SSA_ext)/src/obj/arraySSA/SSA_extension.o $(SSA_ext)/src/obj/utils/sympyComp.o $(SSA_ext)/src/obj/arraySSA/reachingDef_ext.o $(SSA_ext)/src/obj/arraySSA/defsAndUsesTraversal_ext.o -L$(LIBDIR) $(EXTRA_LIBS) -L${OMEGA}/omega_lib/obj/ -lomega -L${OMEGA}/code_gen/obj/ -lcodegen -lm $(CHILL)/parser.tab.o $(CHILL)/lex.yy.o $ $(CHILL)/dep.o $(CHILL)/omegatools.o $(CHILL)/irtools.o $(CHILL)/loop.o $(CHILL)/loop_basic.o $(CHILL)/loop_datacopy.o $(CHILL)/loop_unroll.o $(CHILL)/loop_tile.o $(CHILL)/loop_extra.o $(CHILL)/ir_rose.o
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
	$(CPP) -E -MM -I$(ROSE_INCLUDE) $(EXTRA_CFLGS) $(SRCS) > DEPENDS

include DEPENDS
