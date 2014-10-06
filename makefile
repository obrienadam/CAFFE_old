SHELL = /bin/sh

## Set-up the compiler flags

CC_FLAGS = -O3 -no-prec-div -Wall -fno-alias

CXX_FLAGS = $(CC_FLAGS)

## Add include paths

INCLUDE = -I./$(MATH_DIR)
INCLUDE += -I./$(DATA_STRUCTURES_DIR)
INCLUDE += -I./$(RUN_CONTROL_DIR)
INCLUDE += -I./$(MESH_DIR)
INClUDE += -I./$(FIELD_DIR)

## Directories

MODULES_DIR = Modules
MATH_DIR = src/Math
DATA_STRUCTURES_DIR = src/DataStructures
RUN_CONTROL_DIR = src/RunControl
MESH_DIR = src/Domains/Meshes
FIELD_DIR = src/Fields

## Modules

ADVEC_DIFF = caffeAdvectionDiffusion

MODULES += $(ADVEC_DIFF)

## Object files

ADVEC_DIFF_OBJS = $(ADVEC_DIFF).o

MATH_OBJS += Vector3D.o

RUN_CONTROL_OBJS += RunControl.o
RUN_CONTROL_OBJS += ArgsList.o
RUN_CONTROL_OBJS += Input.o


install: all

all: $(MODULES)

$(ADVEC_DIFF): $(ADVEC_DIFF_OBJS) $(MATH_OBJS) $(RUN_CONTROL_OBJS) HexaFdmMesh.o
	$(CXX) $(INCLUDE) $(CXX_FLAGS) -o $(ADVEC_DIFF) $(ADVEC_DIFF_OBJS) $(MATH_OBJS) $(RUN_CONTROL_OBJS) HexaFdmMesh.o
	mv $(ADVEC_DIFF) bin	

$(ADVEC_DIFF_OBJS):%.o:$(MODULES_DIR)/%.cc
	$(CXX) $(INCLUDE) $(CXX_FLAGS) -c $<

$(MATH_OBJS):%.o:$(MATH_DIR)/%.cc
	$(CXX) $(INCLUDE) $(CXX_FLAGS) -c $<

$(RUN_CONTROL_OBJS):%.o:$(RUN_CONTROL_DIR)/%.cc
	$(CXX) $(INCLUDE) $(CXX_FLAGS) -c $<

HexaFdmMesh.o: $(MESH_DIR)/HexaFdmMesh.cc
	$(CXX) $(INCLUDE) $(CXX_FLAGS) -c $<

clean:
	rm -f *.o *~

superclean: clean
	rm -f bin/*
	rm -f build/*
	rm -f $(MODULES_DIR)/*~
	rm -f $(MATH_DIR)/*~
	rm -f $(DATA_STRUCTURES_DIR)/*~
	rm -f $(RUN_CONTROL_DIR)/*~
	rm -f $(MESH_DIR)/*~
