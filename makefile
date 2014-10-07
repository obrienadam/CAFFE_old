SHELL = /bin/sh

## Set-up the compiler flags

CC_FLAGS = -O3 -no-prec-div -Wall -fno-alias

CXX_FLAGS = $(CC_FLAGS)

## Modules

ADVEC_DIFF = caffeAdvectionDiffusion
MODULES += $(ADVEC_DIFF)

## Directories

BUILD_DIR = Build
MODULES_DIR = Modules
MATH_DIR = src/Math
DATA_STRUCTURES_DIR = src/DataStructures
RUN_CONTROL_DIR = src/RunControl
MESH_DIR = src/Domains/Meshes
FIELD_DIR = src/Fields

ALL_DIRS = $(MODULES_DIR) $(MATH_DIR) $(DATA_STRUCTURES_DIR) $(RUN_CONTROL_DIR) $(MESH_DIR) $(FIELD_DIR)

## Includes

INCLUDE = $(addprefix -I./, $(MATH_DIR) $(DATA_STRUCTURES_DIR) $(RUN_CONTROL_DIR) $(MESH_DIR) $(FIELD_DIR))

## External libraries

## Source files

# Math

MATH_SRC_FILES += Vector3D.cc
MATH_SRC_FILES += Tensor3D.cc

MATH_SRC = $(addprefix $(MATH_DIR)/, $(MATH_SRC_FILES))

# Run control

RUN_CONTROL_SRC_FILES += RunControl.cc
RUN_CONTROL_SRC_FILES += ArgsList.cc
RUN_CONTROL_SRC_FILES += Input.cc

RUN_CONTROL_SRC = $(addprefix $(RUN_CONTROL_DIR)/, $(RUN_CONTROL_SRC_FILES))

# Meshes

MESH_SRC_FILES += PrimitiveMesh.cc
MESH_SRC_FILES += HexaFdmMesh.cc

MESH_SRC = $(addprefix $(MESH_DIR)/, $(MESH_SRC_FILES))

## Module Dependencies

# Advection Diffusion

ADVEC_DIFF_SRC += Modules/$(ADVEC_DIFF).cc
ADVEC_DIFF_SRC += $(MATH_SRC)
ADVEC_DIFF_SRC += $(RUN_CONTROL_SRC)
ADVEC_DIFF_SRC += $(MESH_SRC)
ADVEC_DIFF_OBJS = $(ADVEC_DIFF_SRC:.cc=.o)

install: all

all: $(MODULES)

$(ADVEC_DIFF): $(ADVEC_DIFF_OBJS)
	$(CXX) $(INCLUDE) $(CXX_FLAGS) -o $(ADVEC_DIFF) $(ADVEC_DIFF_OBJS) -lboost_program_options
	mv $(ADVEC_DIFF) bin/

$(ADVEC_DIFF_OBJS):%.o: %.cc
	$(CXX) $(INCLUDE) $(CXX_FLAGS) -c $< -o $@

clean:
	rm -f $(addsuffix /*.o, $(ALL_DIRS))
	rm -f $(addsuffix /*~, $(ALL_DIRS))
