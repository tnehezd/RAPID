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
# Generate include flags for subdirectories
INC_FLAGS = -I$(INC_DIR) $(addprefix -I$(INC_DIR)/, $(SUBDIRS))

LIBOMP := /usr/local/opt/libomp

CFLAGS = -Wall -Wextra -std=c99 -g -O0 $(INC_FLAGS) -D_GNU_SOURCE \
         -I$(LIBOMP)/include -Xpreprocessor -fopenmp

CXXFLAGS = -Wall -Wextra -std=c++17 -g -O0 $(INC_FLAGS) \
           -I$(EXTERN_DIR)/../include -D_GNU_SOURCE \
           -I$(LIBOMP)/include -Xpreprocessor -fopenmp

LDFLAGS = -lm -L$(LIBOMP)/lib -lomp \
          -L/usr/local/Cellar/hdf5/2.0.0_1/lib -lhdf5 -lhdf5_hl

# Recursively find all .c and .cpp source files in the src and extern/src directories
SRCS_C = $(shell find $(SRC_DIR) -name "*.c")
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