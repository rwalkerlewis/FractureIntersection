#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <set>
#include "fracture_element.h"
#include "spatial_grid.h"
#include "intersection.h"

std::vector<FractureElement> generateFractureElements(
    const std::vector<double>& centers, 
    const std::vector<double>& orientations_deg, 
    double length, 
    double height) {
    
        std::vector<FractureElement> fractures;
        
        for (size_t i = 0; i < centers.size() / 3; ++i) {
            Point<3> center(centers[3*i], centers[3*i + 1], centers[3*i + 2]);
            double theta = orientations_deg[i];
            fractures.emplace_back(center, theta, length, height, static_cast<int>(i));
        }
    return fractures;
}

double computeOptimalCellSize(const std::vector<FractureElement>& fractures) {
    double max_size = 0.0;
    for (const auto& fracture : fractures) {
        max_size = std::max(max_size, std::max(fracture.getHeight(), fracture.getLength()));
    }
    return max_size;
}

// Function to draw ASCII map of fractures in XY plane
void drawASCIIMap(const std::vector<FractureElement>& elements, 
                  double map_width = 30.0, double map_height = 20.0, 
                  int char_width = 80, int char_height = 30) {
    
    // Find bounds
    double min_x = -map_width / 2.0, max_x = map_width / 2.0;
    double min_y = -map_height / 2.0, max_y = map_height / 2.0;
    
    // Create ASCII grid
    std::vector<std::string> grid(char_height, std::string(char_width, ' '));
    
    // Draw each fracture
    for (size_t f = 0; f < elements.size(); ++f) {
        const auto& elem = elements[f];
        Point<2> p1, p2;
        elem.getHorizontalEndpoints(p1, p2);
        
        // Convert to grid coordinates
        auto toGridX = [&](double x) -> int {
            return static_cast<int>((x - min_x) / (max_x - min_x) * (char_width - 1));
        };
        auto toGridY = [&](double y) -> int {
            return char_height - 1 - static_cast<int>((y - min_y) / (max_y - min_y) * (char_height - 1));
        };
        
        int x1 = toGridX(p1.x());
        int y1 = toGridY(p1.y());
        int x2 = toGridX(p2.x());
        int y2 = toGridY(p2.y());
        
        // Clamp to grid bounds
        x1 = std::max(0, std::min(char_width - 1, x1));
        y1 = std::max(0, std::min(char_height - 1, y1));
        x2 = std::max(0, std::min(char_width - 1, x2));
        y2 = std::max(0, std::min(char_height - 1, y2));
        
        // Bresenham-like line drawing
        int dx = std::abs(x2 - x1);
        int dy = std::abs(y2 - y1);
        int sx = (x2 > x1) ? 1 : -1;
        int sy = (y2 > y1) ? 1 : -1;
        int err = dx - dy;
        
        int x = x1, y = y1;
        char fracture_char = '0' + (f % 10);
        
        while (true) {
            if (x >= 0 && x < char_width && y >= 0 && y < char_height) {
                if (grid[y][x] == ' ') {
                    grid[y][x] = fracture_char;
                } else if (grid[y][x] != fracture_char) {
                    grid[y][x] = 'X';  // Mark intersections
                }
            }
            
            if (x == x2 && y == y2) break;
            
            int e2 = 2 * err;
            if (e2 > -dy) {
                err -= dy;
                x += sx;
            }
            if (e2 < dx) {
                err += dx;
                y += sy;
            }
        }
    }
    
    // Draw border and axes
    for (int y = 0; y < char_height; ++y) {
        grid[y][0] = (y == char_height / 2) ? '-' : '|';
        grid[y][char_width - 1] = (y == char_height / 2) ? '-' : '|';
    }
    for (int x = 0; x < char_width; ++x) {
        if (grid[0][x] == ' ') grid[0][x] = '_';
        if (grid[char_height - 1][x] == ' ') grid[char_height - 1][x] = '_';
        if (grid[char_height / 2][x] == ' ') grid[char_height / 2][x] = '-';
        grid[char_height / 2][0] = '+';
        grid[char_height / 2][char_width - 1] = '+';
    }
    
    // Print grid
    std::cout << "ASCII Map (XY Plane):\n";
    for (const auto& row : grid) {
        std::cout << row << "\n";
    }
    
    // Print legend
    std::cout << "\nLegend:\n";
    for (size_t i = 0; i < std::min(elements.size(), (size_t)10); ++i) {
        std::cout << "  " << static_cast<char>('0' + i) << " = Fracture " << i << " (θ=" 
                  << elements[i].getTheta() << "°)\n";
    }
    std::cout << "  X = Intersection point\n";
    std::cout << "  + = Origin\n";
}

// Function to draw ASCII map of the spatial grid in XY plane (summed across all Z)
void drawSpatialGridMap(const SpatialGrid& grid, double cell_size, 
                        int grid_width = 60, int grid_height = 20) {
    
    int min_x, max_x, min_y, max_y, min_z, max_z;
    grid.getGridBounds(min_x, max_x, min_y, max_y, min_z, max_z);
    
    if (min_x > max_x) {
        std::cout << "Grid is empty.\n";
        return;
    }
    
    int grid_x_range = max_x - min_x + 1;
    int grid_y_range = max_y - min_y + 1;
    int grid_z_range = max_z - min_z + 1;
    
    // Create a grid matching the actual number of cells
    std::vector<std::string> ascii_grid(grid_y_range, std::string(grid_x_range * 2, ' '));
    
    // Fill in the grid data
    for (int cell_y = min_y; cell_y <= max_y; ++cell_y) {
        for (int cell_x = min_x; cell_x <= max_x; ++cell_x) {
            std::set<int> all_fractures;
            
            // Aggregate fractures across all Z levels for this (X,Y) cell
            for (int cell_z = min_z; cell_z <= max_z; ++cell_z) {
                std::vector<int> fractures = grid.getFracturesInCell(cell_x, cell_y, cell_z);
                for (int fid : fractures) {
                    all_fractures.insert(fid);
                }
            }
            
            // Map to display coordinates
            int disp_y = max_y - cell_y;  // Flip Y for top-down view
            int disp_x = (cell_x - min_x) * 2;
            
            if (disp_x >= 0 && disp_x + 1 < static_cast<int>(ascii_grid[disp_y].size()) && 
                disp_y >= 0 && disp_y < static_cast<int>(ascii_grid.size())) {
                
                if (all_fractures.empty()) {
                    ascii_grid[disp_y][disp_x] = '.';
                    ascii_grid[disp_y][disp_x + 1] = '.';
                } else if (all_fractures.size() == 1) {
                    int frac_id = *all_fractures.begin();
                    char ch = '0' + (frac_id % 10);
                    ascii_grid[disp_y][disp_x] = ch;
                    ascii_grid[disp_y][disp_x + 1] = ch;
                } else {
                    ascii_grid[disp_y][disp_x] = '[';
                    char ch = '*';
                    if (all_fractures.size() <= 9) {
                        ch = '0' + all_fractures.size();
                    }
                    ascii_grid[disp_y][disp_x + 1] = ch;
                }
            }
        }
    }
    
    // Print results
    std::cout << "Spatial Grid Map (XY Plane, aggregated across all Z levels):\n";
    std::cout << "Cell size: " << cell_size << " m\n";
    std::cout << "Grid dimensions: " << grid_x_range << " × " << grid_y_range << " × " << grid_z_range 
              << " cells (X × Y × Z)\n";
    std::cout << "X range: [" << min_x << ", " << max_x << "]\n";
    std::cout << "Y range: [" << min_y << ", " << max_y << "]\n";
    std::cout << "Z range: [" << min_z << ", " << max_z << "]\n";
    
    // Draw grid with cell coordinates
    std::cout << "\nGrid layout (each cell shown as 'XX'):\n";
    
    // Draw top border with X coordinates
    std::cout << "    ";
    for (int x = min_x; x <= max_x; ++x) {
        std::cout << std::setw(3) << x;
    }
    std::cout << "\n";
    
    // Draw grid with Y coordinates
    for (int y = max_y; y >= min_y; --y) {
        std::cout << std::setw(2) << y << " |";
        for (const char ch : ascii_grid[max_y - y]) {
            std::cout << ch;
        }
        std::cout << "|\n";
    }
    
    // Draw bottom border
    std::cout << "    ";
    for (int x = min_x; x <= max_x; ++x) {
        std::cout << std::setw(3) << x;
    }
    std::cout << "\n";
    
    std::cout << "\nLegend:\n";
    for (size_t i = 0; i < 10; ++i) {
        std::cout << "  " << static_cast<char>('0' + i) << " = Fracture " << i << " (or count in cell)\n";
    }
    std::cout << "  .. = Empty cell\n";
    std::cout << "  [N = Multiple fractures (N = fracture count)\n";
}

int main() {

    // Example Data
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "EXAMPLE: Multiple intersecting fractures\n";
    std::cout << std::string(60, '=') << "\n";
    
    // Four fractures positioned to intersect each other
    // All at the same Z level but different XY positions to ensure intersection
    std::vector<double> centers = {
        -5.0,  0.0,  0.0,   // F0: left side, vertical
        5.0,   0.0,  0.0,   // F1: right side, vertical
        0.0,  -5.0,  0.0,   // F2: bottom, horizontal
        0.0,   5.0,  0.0    // F3: top, horizontal
    };
    
    std::vector<double> orientations_deg = {
        90.0,    // F0: vertical orientation
        90.0,    // F1: vertical orientation
        0.0,     // F2: horizontal orientation
        0.0      // F3: horizontal orientation
    };
    
    double length = 12.0;
    double height = 6.0;
    
    std::cout << "Configuration:\n";
    std::cout << "  4 fractures designed to intersect:\n";
    std::cout << "    F0: center=(-5, 0, 0), θ=90°\n";
    std::cout << "    F1: center=(5, 0, 0), θ=90°\n";
    std::cout << "    F2: center=(0, -5, 0), θ=0°\n";
    std::cout << "    F3: center=(0, 5, 0), θ=0°\n";
    std::cout << "  L=" << length << ", H=" << height << "\n";

    // Validate input data
    if (centers.size() % 3 != 0) {
        std::cerr << "Error: Centers vector size must be a multiple of 3 (x, y, z for each fracture).\n";
        return 1;
    }
    
    if (orientations_deg.size() != centers.size() / 3) {
        std::cerr << "Error: Orientations vector size must match the number of fractures (centers.size() / 3).\n";
        return 1;
    }

    // Generate fracture elements
    auto elements = generateFractureElements(centers, orientations_deg, length, height);

    double cell_size = computeOptimalCellSize(elements);

    size_t num_fractures = elements.size();

    // Initialize connections list
    std::vector<std::vector<int>> connections(num_fractures);

    // Build spatial grid
    SpatialGrid grid(cell_size);

    for (size_t i = 0; i < elements.size(); ++i) {
        grid.insert(elements[i]);
    }

    for (size_t i = 0; i < elements.size(); ++i) {
        const auto& element_i = elements[i];

        // Query spatial grid for nearby fractures
        auto nearby = grid.queryNearby(element_i);

        // Test nearby fractures for potential intersection
        for (const auto* element_j : nearby) {
            int j = element_j->getID();

            // Only process pairs (i, j) where i < j to avoid duplicates
            if (i >= static_cast<size_t>(j)) {
                continue;
            }

            // Broad phase: Check if bounding boxes intersect
            if (!element_i.getBoundingBox().intersects(element_j->getBoundingBox())) {
                continue;  // Skip if bounding boxes don't intersect
            }

            // Narrow phase: Check if the two fractures actually intersect
            if (Intersection::fracturesIntersect(element_i, *element_j)) {
                connections[i].push_back(j);
                connections[j].push_back(i);   
            }
        }
    }




    // Print all elements
    std::cout << "\nGenerated Fracture Elements:\n";
    std::cout << std::string(60, '-') << "\n";
    for (const auto& elem : elements) {
        elem.print();
        std::cout << "\n";
    }
    
    std::cout << "\nAlternative: Print using operator<<:\n";
    std::cout << std::string(60, '-') << "\n";
    for (const auto& elem : elements) {
        std::cout << elem << "\n";
    }
    
    // Print connections
    std::cout << "\nFracture Connections (Intersections):\n";
    std::cout << std::string(60, '-') << "\n";
    bool has_connections = false;
    for (size_t i = 0; i < connections.size(); ++i) {
        if (!connections[i].empty()) {
            has_connections = true;
            std::cout << "Fracture " << i << " intersects with: ";
            for (size_t j = 0; j < connections[i].size(); ++j) {
                std::cout << connections[i][j];
                if (j < connections[i].size() - 1) {
                    std::cout << ", ";
                }
            }
            std::cout << "\n";
        }
    }
    if (!has_connections) {
        std::cout << "No fracture intersections detected.\n";
    }
    
    // Print spatial grid statistics
    std::cout << "\nSpatial Grid Statistics:\n";
    std::cout << std::string(60, '-') << "\n";
    int num_cells, num_fractures_in_grid;
    double fill_ratio;
    grid.getGridStatistics(num_cells, num_fractures_in_grid, fill_ratio);
    
    Point<3> domain_min, domain_max;
    grid.getDomainBounds(domain_min, domain_max);
    
    double domain_width = domain_max.x() - domain_min.x();
    double domain_height = domain_max.y() - domain_min.y();
    double domain_depth = domain_max.z() - domain_min.z();
    
    std::cout << "Domain bounds:\n";
    std::cout << "  X: [" << domain_min.x() << ", " << domain_max.x() << "] (width: " << domain_width << ")\n";
    std::cout << "  Y: [" << domain_min.y() << ", " << domain_max.y() << "] (height: " << domain_height << ")\n";
    std::cout << "  Z: [" << domain_min.z() << ", " << domain_max.z() << "] (depth: " << domain_depth << ")\n";
    std::cout << "Cell size: " << cell_size << "\n";
    std::cout << "Occupied grid cells: " << num_cells << "\n";
    std::cout << "Unique fractures: " << num_fractures_in_grid << "\n";
    std::cout << "Fill ratio: " << std::fixed << std::setprecision(2) << (fill_ratio * 100.0) << "% "
              << "(ratio of occupied cells to theoretical grid)\n";
    
    // Draw ASCII map of fractures using actual domain bounds
    std::cout << "\n";
    double map_width = domain_width * 1.1;  // Add 10% padding
    double map_height = domain_height * 1.1;
    drawASCIIMap(elements, map_width, map_height, 100, 40);
    
    // Draw ASCII map of spatial grid
    std::cout << "\n";
    drawSpatialGridMap(grid, cell_size, 50, 40);
    
    return 0;
}
