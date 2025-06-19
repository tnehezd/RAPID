# Compiler to use
CC = gcc

# Compiler flags:
# -O2: Optimization level 2
# -Wall: Enable all warnings
# -Wextra: Enable extra warnings
# -std=c11: Use C11 standard
# -g: Include debugging information
# -lm: Link with the math library
CFLAGS = -O2 -Wall -Wextra -std=c11 -g -lm

# Source files (all .c files in the src directory)
SRCS = $(wildcard src/*.c)

# Object files directory
OBJ_DIR = obj

# Object files (compiled .o files)
OBJS = $(patsubst src/%.c, $(OBJ_DIR)/%.o, $(SRCS))

# Executable name and directory
BIN_DIR = bin
TARGET = $(BIN_DIR)/simulation

# Default target: builds the simulation executable
all: $(TARGET)

# Rule to build the executable
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) -o $@
	@echo "Build successful! Executable in $(TARGET)"

# Rule to compile each .c file into an .o file
$(OBJ_DIR)/%.o: src/%.c
	@mkdir -p $(OBJ_DIR) # Ensure the obj directory exists
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up compiled files and directories
clean:
	@rm -rf $(OBJ_DIR) $(BIN_DIR)
	@echo "Cleaned up build artifacts."

# Phony targets to prevent conflicts with files of the same name
.PHONY: all clean
