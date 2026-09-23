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
OBJ_DIR = obj

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
    
    CFLAGS = -Wall -Wextra -std=c99 -g -O0 $(INC_FLAGS) -D_GNU_SOURCE \
             -I$(LIBOMP)/include -I$(HDF5_HOME)/include -Xpreprocessor -fopenmp
             
    CXXFLAGS = -Wall -Wextra -std=c++17 -g -O0 $(INC_FLAGS) \
               -I$(LIBOMP)/include -I$(HDF5_HOME)/include -Xpreprocessor -fopenmp
               
    LDFLAGS = -lm -L$(LIBOMP)/lib -lomp \
              -L$(HDF5_HOME)/lib -lhdf5 -lhdf5_hl

else
    # Linux (Ubuntu) settings - use gcc/g++ for native OpenMP (libgomp) support
    CC = gcc
    CXX = g++
    
    CFLAGS = -Wall -Wextra -std=c99 -g -O0 $(INC_FLAGS) -D_GNU_SOURCE \
             -I/usr/include/hdf5/serial -fopenmp
             
    CXXFLAGS = -Wall -Wextra -std=c++17 -g -O0 $(INC_FLAGS) \
               -I/usr/include/hdf5/serial -fopenmp
               
    LDFLAGS = -lm -fopenmp \
              -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_hl
endif

# Recursively find all .c files in the src direcotry
# Find all .c files in src, but EXCLUDE python_interface.c
SRCS_C = $(shell find $(SRC_DIR) -name "*.c" ! -name "python_interface.c")

# Generate object file names based on source files
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS_C)) 

# Dependency files (.d)
DEPS = $(OBJS:.o=.d)
DEPS := $(DEPS:.opp=.d)

.PHONY: all clean run debug

all: $(BIN_DIR)/simulation

# Linker
$(BIN_DIR)/simulation: $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@

# C translation (from src folder)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@


-include $(DEPS)

clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR)

run: all
	./$(BIN_DIR)/simulation

debug: all
	lldb ./$(BIN_DIR)/simulation