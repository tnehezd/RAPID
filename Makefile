CC = clang
CXX = clang++

# Könyvtárak
INC_DIR = include
EXTERN_DIR = extern/src
SRC_DIR = src
BIN_DIR = bin
OBJ_DIR = obj

# --- INCLUDE FLAGS ---
SUBDIRS = io physics utils core
INC_FLAGS = -I$(INC_DIR) $(addprefix -I$(INC_DIR)/, $(SUBDIRS)) -I$(INC_DIR)/test

LIBOMP := /usr/local/opt/libomp

CFLAGS = -Wall -Wextra -std=c99 -g -O0 $(INC_FLAGS) -D_GNU_SOURCE \
         -I$(LIBOMP)/include -Xpreprocessor -fopenmp

CXXFLAGS = -Wall -Wextra -std=c++17 -g -O0 $(INC_FLAGS) \
           -I$(EXTERN_DIR)/../include -D_GNU_SOURCE \
           -I$(LIBOMP)/include -Xpreprocessor -fopenmp

LDFLAGS = -lm -L$(LIBOMP)/lib -lomp \
          -L/usr/local/Cellar/hdf5/2.0.0_1/lib -lhdf5 -lhdf5_hl

# --- SOURCE FILE DISCOVERY ---
SRCS_C = $(shell find $(SRC_DIR) -name "*.c")
SRCS_CAPITALC = $(shell find $(SRC_DIR) -name "*.C")
SRCS_CC       = $(shell find $(SRC_DIR) -name "*.cc")
SRCS_CPP      = $(shell find $(SRC_DIR) -name "*.cpp")
EXTERN_CPP    = $(wildcard $(EXTERN_DIR)/*.cpp)

# --- OBJECT FILES ---
OBJS = \
  $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS_C)) \
  $(patsubst $(SRC_DIR)/%.C, $(OBJ_DIR)/%.o, $(SRCS_CAPITALC)) \
  $(patsubst $(SRC_DIR)/%.cc, $(OBJ_DIR)/%.o, $(SRCS_CC)) \
  $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.opp, $(SRCS_CPP)) \
  $(patsubst $(EXTERN_DIR)/%.cpp, $(OBJ_DIR)/%.opp, $(EXTERN_CPP))


# --- DEPENDENCIES ---
DEPS = $(OBJS:.o=.d)
DEPS := $(DEPS:.opp=.d)

.PHONY: all clean run debug test

# --- MAIN SIMULATION BUILD ---
all: $(BIN_DIR)/simulation

$(BIN_DIR)/simulation: $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(OBJS) $(LDFLAGS) -o $@

# --- COMPILATION RULES ---
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.C
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.cc
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

$(OBJ_DIR)/%.opp: $(SRC_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

$(OBJ_DIR)/%.opp: $(EXTERN_DIR)/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# --- CLEAN ---
clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR)

# --- RUN ---
run: all
	./$(BIN_DIR)/simulation

debug: all
	lldb ./$(BIN_DIR)/simulation

-include $(DEPS)
