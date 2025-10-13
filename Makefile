CC = gcc
CFLAGS = -Wall -Wextra -std=c99 -g -O0 -fopenmp -I$(INC_DIR) -D_GNU_SOURCE
LDFLAGS = -lm -fopenmp -g

BIN_DIR = bin
OBJ_DIR = obj
SRC_DIR = src
INC_DIR = include

SRCS = $(wildcard $(SRC_DIR)/*.c)
OBJS = $(patsubst $(SRC_DIR)/%.c, $(OBJ_DIR)/%.o, $(SRCS))

.PHONY: all clean run debug

all: $(BIN_DIR)/simulation

$(BIN_DIR)/simulation: $(OBJS)
	@mkdir -p $(BIN_DIR)
	$(CC) $(CFLAGS) $(OBJS) $(LDFLAGS) -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.c
	@mkdir -p $(OBJ_DIR)
	$(CC) $(CFLAGS) -MMD -MP -c $< -o $@

-include $(OBJS:.o=.d)

clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR)

run: all
	./$(BIN_DIR)/simulation

debug: all
	lldb ./$(BIN_DIR)/simulation
