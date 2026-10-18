# --- PLATFORM-SPECIFIC COMPILER CONFIGURATION ---
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Darwin)
    CC = clang
    CXX = clang++
else
    CC = gcc
    CXX = g++
endif

# Directories
INC_DIR = include
SRC_DIR = src
BIN_DIR = bin
OBJ_ROOT_DIR = obj
BUILD ?= release
.DEFAULT_GOAL := all

ifeq ($(BUILD),debug)
    OPT_FLAGS = -O0 -g
else ifeq ($(BUILD),release)
    OPT_FLAGS = -O2 -g
else
    $(error Unknown BUILD '$(BUILD)'; use debug or release)
endif

OBJ_DIR = $(OBJ_ROOT_DIR)/$(BUILD)
BUILD_BIN_DIR = $(BIN_DIR)/$(BUILD)

# Subdirectories for include files
SUBDIRS = io physics utils core benchmark_tests
INC_FLAGS = -I$(INC_DIR) $(addprefix -I$(INC_DIR)/, $(SUBDIRS))

# --- PLATFORM-SPECIFIC CONFIGURATION ---
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Darwin)
    # Automatically detect Homebrew path (handles both native Apple Silicon and migrated /usr/local setups)
    ifneq ($(wildcard /opt/homebrew/opt/libomp),)
        LIBOMP = /opt/homebrew/opt/libomp
    else
        LIBOMP = /usr/local/opt/libomp
    endif

    ifneq ($(wildcard /opt/homebrew/opt/hdf5),)
        HDF5_HOME = /opt/homebrew/opt/hdf5
    else
        HDF5_HOME = /usr/local/opt/hdf5
    endif
    
    CFLAGS = -Wall -Wextra -std=c99 $(OPT_FLAGS) $(INC_FLAGS) -D_GNU_SOURCE \
             -I$(LIBOMP)/include -I$(HDF5_HOME)/include -Xpreprocessor -fopenmp
             
    CXXFLAGS = -Wall -Wextra -std=c++17 $(OPT_FLAGS) $(INC_FLAGS) \
               -I$(LIBOMP)/include -I$(HDF5_HOME)/include -Xpreprocessor -fopenmp
               
    LDFLAGS = -lm -L$(LIBOMP)/lib -lomp \
              -L$(HDF5_HOME)/lib -lhdf5 -lhdf5_hl

else
    CC = gcc
    CXX = g++

    # Base flags
    CFLAGS = -Wall -Wextra -std=c99 $(OPT_FLAGS) $(INC_FLAGS) -D_GNU_SOURCE -fopenmp
    CXXFLAGS = -Wall -Wextra -std=c++17 $(OPT_FLAGS) $(INC_FLAGS) -fopenmp

    # HDF5 discovery: user override > pkg-config > Ubuntu serial paths
    ifndef HDF5_INCLUDE
    ifndef HDF5_LIB
        HDF5_PC_CFLAGS := $(shell pkg-config --cflags hdf5 2>/dev/null || pkg-config --cflags hdf5-serial 2>/dev/null)
        HDF5_PC_LIBS := $(shell pkg-config --libs-only-L hdf5 2>/dev/null || pkg-config --libs-only-L hdf5-serial 2>/dev/null)
        ifneq ($(strip $(HDF5_PC_CFLAGS)),)
            CFLAGS += $(HDF5_PC_CFLAGS)
            CXXFLAGS += $(HDF5_PC_CFLAGS)
            LDFLAGS += $(HDF5_PC_LIBS)
        else
            ifneq ($(wildcard /usr/include/hdf5/serial/hdf5.h),)
                HDF5_INCLUDE = /usr/include/hdf5/serial
                HDF5_LIB = /usr/lib/$(shell $(CC) -dumpmachine)/hdf5/serial
            endif
        endif
    endif
    endif

    ifdef HDF5_INCLUDE
        CFLAGS += -I$(HDF5_INCLUDE)
        CXXFLAGS += -I$(HDF5_INCLUDE)
    endif

    ifdef HDF5_LIB
        LDFLAGS += -L$(HDF5_LIB)
    endif

    # Always link HDF5
    LDFLAGS += -lm -fopenmp -lhdf5 -lhdf5_hl
endif


# Recursively find all .c files in the src direcotry
# Find all .c files in src, but EXCLUDE python_interface.c
SRCS_C = $(shell find $(SRC_DIR) -name "*.c" ! -name "python_interface.c")

# Generate object file names based on source files
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS_C)) 

# Dependency files (.d)
DEPS = $(OBJS:.o=.d)
DEPS := $(DEPS:.opp=.d)

.PHONY: all clean run debug debug-build release test-thermal test-boundary FORCE

all: $(BIN_DIR)/simulation

# Keep the wrapper's executable path while allowing safe profile switches.
$(BIN_DIR)/simulation: $(BUILD_BIN_DIR)/simulation FORCE
	@cmp -s $< $@ || cp $< $@

# Linker
$(BUILD_BIN_DIR)/simulation: $(OBJS)
	@mkdir -p $(BUILD_BIN_DIR)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@

# C translation (from src folder)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@


-include $(DEPS)

clean:
	@rm -rf $(OBJ_ROOT_DIR) $(BIN_DIR)

run: all
	./$(BIN_DIR)/simulation

release:
	$(MAKE) BUILD=release all

debug-build:
	$(MAKE) BUILD=debug all

debug: debug-build
	lldb ./$(BIN_DIR)/debug/simulation

# Standalone thermal tests need neither HDF5 nor OpenMP.
$(BUILD_BIN_DIR)/test_thermal_model: tests/test_thermal_model.c src/physics/thermal_model.c src/physics/irradiation.c $(wildcard include/physics/*.h include/core/*.h)
	@mkdir -p $(BUILD_BIN_DIR)
	$(CC) -Wall -Wextra -Werror -std=c99 $(OPT_FLAGS) $(INC_FLAGS) \
		tests/test_thermal_model.c src/physics/thermal_model.c src/physics/irradiation.c -lm -o $@

$(BIN_DIR)/test_thermal_model: $(BUILD_BIN_DIR)/test_thermal_model FORCE
	@cmp -s $< $@ || cp $< $@

test-thermal: $(BIN_DIR)/test_thermal_model
	./$(BIN_DIR)/test_thermal_model

$(BUILD_BIN_DIR)/test_boundary_conditions: tests/test_boundary_conditions.c src/physics/boundary_conditions.c src/physics/thermal_model.c src/physics/irradiation.c $(wildcard include/physics/*.h include/core/*.h include/io/*.h)
	@mkdir -p $(BUILD_BIN_DIR)
	$(CC) -Wall -Wextra -Werror -std=c99 $(OPT_FLAGS) $(INC_FLAGS) \
		tests/test_boundary_conditions.c src/physics/boundary_conditions.c \
		src/physics/thermal_model.c src/physics/irradiation.c -lm -o $@

test-boundary: $(BUILD_BIN_DIR)/test_boundary_conditions
	./$(BUILD_BIN_DIR)/test_boundary_conditions

$(BUILD_BIN_DIR)/test_energy_equation: tests/test_energy_equation.c $(filter-out $(OBJ_DIR)/core/main.o,$(OBJS))
	@mkdir -p $(BUILD_BIN_DIR)
	$(CC) $(CFLAGS) -c $< -o $(OBJ_DIR)/test_energy_equation.o
	$(CXX) $(OBJ_DIR)/test_energy_equation.o $(filter-out $(OBJ_DIR)/core/main.o,$(OBJS)) $(LDFLAGS) -o $@

.PHONY: test-energy
test-energy: $(BUILD_BIN_DIR)/test_energy_equation
	./$(BUILD_BIN_DIR)/test_energy_equation