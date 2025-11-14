# Makefile for FractureIntersection Project
# Builds the fracture intersection analysis application

# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++17 -Wall -Wextra -Wpedantic
CPPFLAGS := -I./include
LDFLAGS := -lm

# Build flags
DEBUG_FLAGS := -g -O0 -DDEBUG
RELEASE_FLAGS := -O2 -DNDEBUG

# Directories
SRC_DIR := src
INCLUDE_DIR := include
BUILD_DIR := build
BIN_DIR := bin

# Source and object files
SOURCES := $(wildcard $(SRC_DIR)/*.cpp)
OBJECTS := $(SOURCES:$(SRC_DIR)/%.cpp=$(BUILD_DIR)/%.o)
EXECUTABLE := $(BIN_DIR)/fracture_app

# Phony targets
.PHONY: all debug release clean help run test

# Default target
all: release

# Debug build
debug: CXXFLAGS += $(DEBUG_FLAGS)
debug: $(EXECUTABLE)

# Release build
release: CXXFLAGS += $(RELEASE_FLAGS)
release: $(EXECUTABLE)

# Main executable
$(EXECUTABLE): $(OBJECTS) | $(BIN_DIR)
	@echo "Linking $@..."
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(OBJECTS) -o $@ $(LDFLAGS)
	@echo "Build complete: $@"

# Object files
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.cpp | $(BUILD_DIR)
	@echo "Compiling $<..."
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) -c $< -o $@

# Create directories
$(BUILD_DIR) $(BIN_DIR):
	@mkdir -p $@

# Run the application
run: release
	@echo "Running fracture_app..."
	@./$(EXECUTABLE)

# Compile test program
test: $(BIN_DIR)/test_bbox
	@echo "Running bbox test..."
	@./$(BIN_DIR)/test_bbox

$(BIN_DIR)/test_bbox: test_bbox.cpp $(INCLUDE_DIR)/fracture_element.h $(INCLUDE_DIR)/geometry.h | $(BIN_DIR)
	@echo "Compiling test_bbox..."
	@$(CXX) $(CXXFLAGS) $(CPPFLAGS) test_bbox.cpp -o $@

# Clean build artifacts
clean:
	@echo "Cleaning build artifacts..."
	@rm -rf $(BUILD_DIR) $(BIN_DIR)
	@echo "Clean complete"

# Help target
help:
	@echo "FractureIntersection Makefile"
	@echo ""
	@echo "Targets:"
	@echo "  all        Build in release mode (default)"
	@echo "  debug      Build in debug mode with symbols"
	@echo "  release    Build optimized release version"
	@echo "  run        Build and run the application"
	@echo "  test       Build and run bbox test"
	@echo "  clean      Remove all build artifacts"
	@echo "  help       Show this help message"
	@echo ""
	@echo "Examples:"
	@echo "  make              # Build release version"
	@echo "  make debug        # Build debug version"
	@echo "  make run          # Build and run application"
	@echo "  make clean        # Clean build directory"
	@echo "  make test         # Run unit tests"
