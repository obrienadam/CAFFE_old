SHELL = /bin/sh

## Set-up the compiler flags

CC_FLAGS = -O3 -Wall -std=c++11

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
DOMAIN_DIR = src/Domains
MESH_DIR = src/Domains/Meshes
STATE_DIR = src/State
SOLVER_DIR = src/Solvers
SCHEME_DIR = src/Schemes

ALL_DIRS = $(MODULES_DIR) $(MATH_DIR) $(DATA_STRUCTURES_DIR) $(RUN_CONTROL_DIR) \
$(DOMAIN_DIR) $(MESH_DIR) $(STATE_DIR) $(SOLVER_DIR) $(SCHEME_DIR)

## Includes

INCLUDE = $(addprefix -I./, $(ALL_DIRS))

## External libraries

LIBS += -lboost_program_options 
LIBS += -lboost_system 
LIBS += -lboost_date_time 
LIBS += -lboost_chrono

## Source files

# Math

MATH_SRC_FILES += Vector3D.cc
MATH_SRC_FILES += Tensor3D.cc

MATH_SRC = $(addprefix $(MATH_DIR)/, $(MATH_SRC_FILES))

# Run control

RUN_CONTROL_SRC_FILES += RunControl.cc
RUN_CONTROL_SRC_FILES += ArgsList.cc
RUN_CONTROL_SRC_FILES += Input.cc
RUN_CONTROL_SRC_FILES += Output.cc

RUN_CONTROL_SRC = $(addprefix $(RUN_CONTROL_DIR)/, $(RUN_CONTROL_SRC_FILES))

# Schemes

SCHEME_SRC_FILES += FiniteDifference.cc

SCHEME_SRC = $(addprefix $(SCHEME_DIR)/, $(SCHEME_SRC_FILES))

## Module Dependencies

# Advection Diffusion

ADVEC_DIFF_SRC += Modules/$(ADVEC_DIFF).cc
ADVEC_DIFF_SRC += $(STATE_DIR)/AdvectionDiffusion.cc
ADVEC_DIFF_SRC += $(MATH_SRC)
ADVEC_DIFF_SRC += $(RUN_CONTROL_SRC)
ADVEC_DIFF_SRC += $(DOMAIN_SRC)
ADVEC_DIFF_SRC += $(MESH_SRC)
ADVEC_DIFF_SRC += $(SOLVER_SRC)
ADVEC_DIFF_SRC += $(SCHEME_SRC)
ADVEC_DIFF_OBJS = $(ADVEC_DIFF_SRC:.cc=.o)

install: all

all: $(MODULES)

$(ADVEC_DIFF): $(ADVEC_DIFF_OBJS)
	$(CXX) $(INCLUDE) $(CXX_FLAGS) -o $(ADVEC_DIFF) $(ADVEC_DIFF_OBJS) $(LIBS)
	mv $(ADVEC_DIFF) bin/

$(ADVEC_DIFF_OBJS):%.o: %.cc
	$(CXX) $(INCLUDE) $(CXX_FLAGS) -c $< -o $@

clean:
	rm -f $(addsuffix /*.o, $(ALL_DIRS))
	rm -f $(addsuffix /*~, $(ALL_DIRS))
