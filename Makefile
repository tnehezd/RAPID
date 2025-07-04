# Makefile
CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -g -O0 -fopenmp 	# Add -O0 -g for easier debugging
LDFLAGS = -lm -fopenmp -g 		# Added -g here explicitly for linking

BIN_DIR = bin
OBJ_DIR = obj
SRC_DIR = src
INC_DIR = include

# C forrásfájlok listája (mostantól az init_tool is benne van)
SRCS = $(SRC_DIR)/main.c \
       $(SRC_DIR)/config.c \
       $(SRC_DIR)/disk_model.c \
       $(SRC_DIR)/dust_physics.c \
       $(SRC_DIR)/io_utils.c \
       $(SRC_DIR)/simulation_core.c \
       $(SRC_DIR)/utils.c \
       $(SRC_DIR)/init_tool_module.c \
       $(SRC_DIR)/parser.c \
       $(SRC_DIR)/particle_data.c

# Objektumfájlok listája
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

# Header fájlok (csak a függőségekhez)
HEADERS = $(wildcard $(INC_DIR)/*.h) $(SRC_DIR)/init_tool_module.h # Hozzáadjuk az új headert

# Célok
.PHONY: all clean run debug

all: $(BIN_DIR)/simulation

$(BIN_DIR)/simulation: $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS)
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -I$(INC_DIR) -c $< -o $@

clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR)

run: all
	./$(BIN_DIR)/simulation

debug: all
	lldb ./$(BIN_DIR)/simulation