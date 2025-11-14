# FractureIntersection

A C++ application for detecting intersections between planar fracture elements in 3D space using an efficient spatial grid data structure.

## Features

- **3D Fracture Representation**: Planar fractures defined by center point, orientation (theta), length, and height
- **Spatial Grid Acceleration**: Sparse hash-based spatial grid for efficient intersection queries
- **Multiple Test Cases**: 7 example scenarios including random fracture generation
- **Visualization**: ASCII maps of fractures and spatial grid in XY plane
- **Detailed Output**: Timestamped output files with fracture summaries, connection tables, and statistics

## Building

### Requirements
- C++17 compatible compiler (g++, clang++)
- Make

### Compilation

```bash
# Build with default settings (release mode)
make

# Build in debug mode
make debug

# Build in release mode (optimized)
make release

# Clean build artifacts
make clean

# Build and run
make run
```

### Build Targets

- `make` or `make all` - Build the application (release mode by default)
- `make debug` - Build with debug symbols and no optimization
- `make release` - Build with optimizations
- `make run` - Build and execute the application
- `make clean` - Remove all build artifacts
- `make help` - Display available targets

## Running

After building, run the executable:

```bash
./bin/fracture_app
```

Output will be displayed on the console and saved to a timestamped file:
```
fracture_output_YYYYMMDD_HHMMSS.txt
```

## Project Structure

```
FractureIntersection/
├── include/           # Header files
│   ├── fracture_element.h
│   ├── geometry.h
│   ├── intersection.h
│   └── spatial_grid.h
├── src/              # Source files
│   └── main.cpp
├── bin/              # Compiled binaries (generated)
├── build/            # Build artifacts (generated)
├── Makefile
└── README.md
```

## Examples

The application includes 7 test scenarios:

1. **Cross-Pattern Intersections** - Four fractures in a cross pattern
2. **Radial Fracture Pattern** - Six fractures radiating from center
3. **Parallel Fractures** - Non-intersecting parallel fractures
4. **3D Multi-Level Fractures** - Fractures at different Z levels
5. **Dense Cluster** - Ten fractures in close proximity
6. **XY Overlap with Z Separation** - Demonstrates Z-axis separation
7. **Random Fracture Generation** - 30 randomly oriented fractures over a user-defined domain

## Algorithm

The application uses a sparse spatial grid for efficient intersection detection:
- Fractures are inserted into grid cells based on their bounding boxes
- Only occupied cells are stored in memory (hash map based)
- Intersection queries only check fractures in nearby cells
- O(1) cell lookup with significant performance improvement for large datasets
