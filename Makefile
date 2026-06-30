CC = clang
CXX = clang++

INC_DIR = include
EXTERN_DIR = extern/src

BIN_DIR = bin
OBJ_DIR = obj
SRC_DIR = src

LIBOMP := /usr/local/opt/libomp

CFLAGS = -Wall -Wextra -std=c99 -g -O0 \
         -I$(INC_DIR) -D_GNU_SOURCE \
         -I$(LIBOMP)/include \
         -Xpreprocessor -fopenmp

CXXFLAGS = -Wall -Wextra -std=c++17 -g -O0 \
           -I$(INC_DIR) \
           -I$(EXTERN_DIR)/../include \
           -D_GNU_SOURCE \
           -I/usr/local/opt/libomp/include \
           -Xpreprocessor -fopenmp

LDFLAGS = -lm \
          -L$(LIBOMP)/lib -lomp \
          -L/usr/local/Cellar/hdf5/2.0.0_1/lib -lhdf5 -lhdf5_hl

# C forrásfájlok keresése a src/ mappában
SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

# C++ forrásfájlok keresése MINDKÉT mappában (src/ és extern/src/)
SRC_CPP_SRCS = $(wildcard $(SRC_DIR)/*.cpp)
EXTERN_CPP_SRCS = $(wildcard $(EXTERN_DIR)/*.cpp)

# Object fájlok leképezése (.opp kiterjesztéssel)
SRC_CPP_OBJS = $(patsubst $(SRC_DIR)/%.cpp, $(OBJ_DIR)/%.opp, $(SRC_CPP_SRCS))
EXTERN_CPP_OBJS = $(patsubst $(EXTERN_DIR)/%.cpp, $(OBJ_DIR)/%.opp, $(EXTERN_CPP_SRCS))

# Összesített C++ object lista
CPP_OBJS = $(SRC_CPP_OBJS) $(EXTERN_CPP_OBJS)

# Összesített függőségi (.d) fájlok listája az automatikus fejléc-követéshez
DEPS = $(OBJS:.o=.d) $(CPP_OBJS:.opp=.d)

.PHONY: all clean run debug

all: $(BIN_DIR)/simulation

# LINKELÉS
$(BIN_DIR)/simulation: $(OBJS) $(CPP_OBJS)
	@mkdir -p $(BIN_DIR)
	$(CXX) $(OBJS) $(CPP_OBJS) $(LDFLAGS) -o $@

# C FÁJLOK FORDÍTÁSA (src/*.c -> obj/*.o)
$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

# C++ FÁJLOK FORDÍTÁSA A SRC/ MAPPÁBÓL (src/*.cpp -> obj/*.opp)
$(OBJ_DIR)/%.opp: $(SRC_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# C++ FÁJLOK FORDÍTÁSA AZ EXTERN/ MAPPÁBÓL (extern/src/*.cpp -> obj/*.opp)
$(OBJ_DIR)/%.opp: $(EXTERN_DIR)/%.cpp
	@mkdir -p $(OBJ_DIR)
	$(CXX) $(CXXFLAGS) -MMD -MP -c $< -o $@

# Függőségek beágyazása (ha változik egy .h/.hpp, a make tudni fogja, mit kell újrafordítani)
-include $(DEPS)

clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR)

run: all
	./$(BIN_DIR)/simulation

debug: all
	lldb ./$(BIN_DIR)/simulation