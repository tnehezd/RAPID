CC = clang
CXX = clang++

# Könyvtárak
INC_DIR = include
EXTERN_DIR = extern/src
SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = obj

# Subdirectories for include files
SUBDIRS = io physics utils core
INC_FLAGS = -I$(INC_DIR) $(addprefix -I$(INC_DIR)/, $(SUBDIRS))

# --- PLATFORM-FÜGGŐ BEÁLLÍTÁSOK ---
UNAME_S := $(shell uname -s)

ifeq ($(UNAME_S), Darwin)
    # macOS (Homebrew) beállítások
    LIBOMP ?= /usr/local/opt/libomp
    HDF5_HOME ?= /usr/local/opt/hdf5
    
    CFLAGS = -Wall -Wextra -std=c99 -g -O0 $(INC_FLAGS) -D_GNU_SOURCE \
             -I$(LIBOMP)/include -I$(HDF5_HOME)/include -Xpreprocessor -fopenmp
             
    CXXFLAGS = -Wall -Wextra -std=c++17 -g -O0 $(INC_FLAGS) \
               -I$(EXTERN_DIR)/../include -D_GNU_SOURCE \
               -I$(LIBOMP)/include -I$(HDF5_HOME)/include -Xpreprocessor -fopenmp
               
    LDFLAGS = -lm -L$(LIBOMP)/lib -lomp \
              -L$(HDF5_HOME)/lib -lhdf5 -lhdf5_hl
else
    # Linux (Ubuntu / GitHub Actions) beállítások
    CFLAGS = -Wall -Wextra -std=c99 -g -O0 $(INC_FLAGS) -D_GNU_SOURCE \
             -I/usr/include/hdf5/serial -fopenmp
             
    CXXFLAGS = -Wall -Wextra -std=c++17 -g -O0 $(INC_FLAGS) \
               -I$(EXTERN_DIR)/../include -D_GNU_SOURCE \
               -I/usr/include/hdf5/serial -fopenmp
               
    LDFLAGS = -lm -fopenmp \
              -L/usr/lib/x86_64-linux-gnu/hdf5/serial -lhdf5 -lhdf5_hl
endif

# Recursively find all .c and .cpp source files in the src and extern/src directories
# Find all .c files in src, but EXCLUDE python_interface.c
SRCS_C = $(shell find $(SRC_DIR) -name "*.c" ! -name "python_interface.c")
SRCS_CPP = $(shell find $(SRC_DIR) -name "*.cpp")
EXTERN_CPP = $(wildcard $(EXTERN_DIR)/*.cpp)

# Generate object file names based on source files
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS_C)) \
       $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.opp, $(SRCS_CPP)) \
       $(patsubst $(EXTERN_DIR)/%.cpp, $(OBJ_DIR)/%.opp, $(EXTERN_CPP))

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

# C++ translation (from src folder) -- for cpp wrappers and other C++ files
$(OBJ_DIR)/%.opp: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# C++ translation (from extern folder)
$(OBJ_DIR)/%.opp: $(EXTERN_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

-include $(DEPS)

clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR)

run: all
	./$(BIN_DIR)/simulation

debug: all
	lldb ./$(BIN_DIR)/simulation