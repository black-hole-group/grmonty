# Source directory
SRC_DIR = src

# Build directory
BUILD_DIR = build

# Compiler and flags
CC = gcc
CFLAGS = -pg -g -Wall -O0 -fopenmp
LDFLAGS = -lm -lgsl -lgslcblas -fopenmp

# Source files
SRCS = $(SRC_DIR)/grmonty.c $(SRC_DIR)/compton.c $(SRC_DIR)/init_geometry.c \
       $(SRC_DIR)/tetrads.c $(SRC_DIR)/geodesics.c $(SRC_DIR)/radiation.c \
       $(SRC_DIR)/jnu_mixed.c $(SRC_DIR)/hotcross.c \
       $(SRC_DIR)/track_super_photon.c $(SRC_DIR)/scatter_super_photon.c \
       $(SRC_DIR)/harm_model.c $(SRC_DIR)/harm_utils.c \
       $(SRC_DIR)/hamr_model.c $(SRC_DIR)/init_harm_data.c \
       $(SRC_DIR)/init_hamr_data2D.c $(SRC_DIR)/init_hamr_data3D.c 

# Object files
OBJS = $(patsubst $(SRC_DIR)/%.c,$(BUILD_DIR)/%.o,$(SRCS))

# Include files
INCS = $(SRC_DIR)/decs.h $(SRC_DIR)/constants.h $(SRC_DIR)/harm_model.h

# Executable
EXECUTABLE = grmonty

# Build rule
$(EXECUTABLE): $(OBJS) $(INCS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -o $(EXECUTABLE) $(OBJS) $(LDFLAGS)

# Compile rule
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(INCS) | $(BUILD_DIR)
	$(CC) $(CFLAGS) -c -o $@ $<

# Create build directory if it doesn't exist
$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

# Clean rule
clean:
	/bin/rm -rf $(BUILD_DIR) grmonty
