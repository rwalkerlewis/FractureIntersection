#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <vector>
#include <iomanip>
#include <set>
#include <ctime>
#include <random>
#include "fracture_element.h"
#include "spatial_grid.h"
#include "intersection.h"

// Stream buffer that writes to multiple outputs (like Unix 'tee')
class TeeBuf : public std::streambuf {
private:
    std::streambuf* sb1_;
    std::streambuf* sb2_;
    
public:
    TeeBuf(std::streambuf* sb1, std::streambuf* sb2) : sb1_(sb1), sb2_(sb2) {}
    
protected:
    virtual int overflow(int c) override {
        if (c == EOF) {
            return !EOF;
        }
        int const r1 = sb1_->sputc(c);
        int const r2 = sb2_->sputc(c);
        return (r1 == EOF || r2 == EOF) ? EOF : c;
    }
    
    virtual int sync() override {
        int const r1 = sb1_->pubsync();
        int const r2 = sb2_->pubsync();
        return (r1 == 0 && r2 == 0) ? 0 : -1;
    }
};

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

// Function to generate random fractures within a domain
std::vector<FractureElement> generateRandomFractures(
    int num_fractures,
    double domain_x_min, double domain_x_max,
    double domain_y_min, double domain_y_max,
    double domain_z_min, double domain_z_max,
    double min_length, double max_length,
    double min_height, double max_height,
    unsigned int seed = 42) {
    
    std::vector<FractureElement> fractures;
    std::mt19937 gen(seed);
    
    // Uniform distributions for position and properties
    std::uniform_real_distribution<> dist_x(domain_x_min, domain_x_max);
    std::uniform_real_distribution<> dist_y(domain_y_min, domain_y_max);
    std::uniform_real_distribution<> dist_z(domain_z_min, domain_z_max);
    std::uniform_real_distribution<> dist_theta(0.0, 180.0);
    std::uniform_real_distribution<> dist_length(min_length, max_length);
    std::uniform_real_distribution<> dist_height(min_height, max_height);
    
    for (int i = 0; i < num_fractures; ++i) {
        Point<3> center(dist_x(gen), dist_y(gen), dist_z(gen));
        double theta = dist_theta(gen);
        double length = dist_length(gen);
        double height = dist_height(gen);
        
        fractures.emplace_back(center, theta, length, height, i);
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
    
    // Find actual bounds from fracture endpoints
    double min_x = 0.0, max_x = 0.0, min_y = 0.0, max_y = 0.0;
    bool first = true;
    
    for (const auto& elem : elements) {
        Point<2> p1, p2;
        elem.getHorizontalEndpoints(p1, p2);
        
        if (first) {
            min_x = std::min(p1.x(), p2.x());
            max_x = std::max(p1.x(), p2.x());
            min_y = std::min(p1.y(), p2.y());
            max_y = std::max(p1.y(), p2.y());
            first = false;
        } else {
            min_x = std::min(min_x, std::min(p1.x(), p2.x()));
            max_x = std::max(max_x, std::max(p1.x(), p2.x()));
            min_y = std::min(min_y, std::min(p1.y(), p2.y()));
            max_y = std::max(max_y, std::max(p1.y(), p2.y()));
        }
    }
    
    // Add padding (10% on each side)
    double x_range = max_x - min_x;
    double y_range = max_y - min_y;
    double x_padding = x_range * 0.1;
    double y_padding = y_range * 0.1;
    
    min_x -= x_padding;
    max_x += x_padding;
    min_y -= y_padding;
    max_y += y_padding;
    
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
                        int grid_width [[maybe_unused]] = 60, int grid_height [[maybe_unused]] = 20) {
    
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

// Function to run a fracture intersection example
void runExample(const std::string& example_name, 
                const std::string& description,
                const std::vector<double>& centers,
                const std::vector<double>& orientations_deg,
                double length,
                double height) {
    
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "EXAMPLE: " << example_name << "\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << description << "\n\n";
    
    std::cout << "Configuration:\n";
    std::cout << "  Number of fractures: " << (centers.size() / 3) << "\n";
    std::cout << "  Length: " << length << " m, Height: " << height << " m\n\n";

    // Validate input data
    if (centers.size() % 3 != 0) {
        std::cerr << "Error: Centers vector size must be a multiple of 3 (x, y, z for each fracture).\n";
        return;
    }
    
    if (orientations_deg.size() != centers.size() / 3) {
        std::cerr << "Error: Orientations vector size must match the number of fractures (centers.size() / 3).\n";
        return;
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

    // Find intersections
    for (size_t i = 0; i < elements.size(); ++i) {
        const auto& element_i = elements[i];
        auto nearby = grid.queryNearby(element_i);

        for (const auto* element_j : nearby) {
            int j = element_j->getID();
            if (i >= static_cast<size_t>(j)) {
                continue;
            }

            if (!element_i.getBoundingBox().intersects(element_j->getBoundingBox())) {
                continue;
            }

            if (Intersection::fracturesIntersect(element_i, *element_j)) {
                connections[i].push_back(j);
                connections[j].push_back(i);   
            }
        }
    }

    // Print fracture summary with centerpoints
    std::cout << "\nFracture Summary:\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::left << std::setw(4) << "ID" 
              << std::setw(25) << "Center (x, y, z)" 
              << std::setw(10) << "Theta(°)" 
              << std::setw(10) << "Length" 
              << std::setw(10) << "Height" << "\n";
    std::cout << std::string(80, '-') << "\n";
    for (const auto& elem : elements) {
        auto center = elem.getCenter();
        std::cout << std::left << std::setw(4) << elem.getID()
                  << "(" << std::fixed << std::setprecision(2) << std::setw(6) << center.x() 
                  << ", " << std::setw(6) << center.y() 
                  << ", " << std::setw(6) << center.z() << ")  "
                  << std::setw(10) << elem.getTheta()
                  << std::setw(10) << elem.getLength()
                  << std::setw(10) << elem.getHeight() << "\n";
    }
    std::cout << std::string(80, '-') << "\n";

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
    double map_width = domain_width * 1.1;
    double map_height = domain_height * 1.1;
    drawASCIIMap(elements, map_width, map_height, 100, 40);
    
    // Draw ASCII map of spatial grid
    std::cout << "\n";
    drawSpatialGridMap(grid, cell_size, 50, 40);
}

// Overloaded runExample for pre-generated fractures
void runExample(const std::string& example_name,
                const std::string& description,
                const std::vector<FractureElement>& elements) {
    
    std::cout << "\n" << std::string(80, '=') << "\n";
    std::cout << "EXAMPLE: " << example_name << "\n";
    std::cout << std::string(80, '=') << "\n";
    std::cout << description << "\n\n";
    
    std::cout << "Configuration:\n";
    std::cout << "  Number of fractures: " << elements.size() << "\n\n";

    // Print fracture summary table
    std::cout << "Fracture Summary:\n";
    std::cout << std::string(80, '-') << "\n";
    std::cout << std::setw(4) << "ID" << " | "
              << std::setw(10) << "Center X" << " | "
              << std::setw(10) << "Center Y" << " | "
              << std::setw(10) << "Center Z" << " | "
              << std::setw(8) << "Theta" << " | "
              << std::setw(8) << "Length" << " | "
              << std::setw(8) << "Height" << "\n";
    std::cout << std::string(80, '-') << "\n";
    
    for (const auto& frac : elements) {
        std::cout << std::setw(4) << frac.getID() << " | "
                  << std::setw(10) << std::fixed << std::setprecision(2) << frac.getCenter().x() << " | "
                  << std::setw(10) << std::fixed << std::setprecision(2) << frac.getCenter().y() << " | "
                  << std::setw(10) << std::fixed << std::setprecision(2) << frac.getCenter().z() << " | "
                  << std::setw(8) << std::fixed << std::setprecision(1) << frac.getTheta() << " | "
                  << std::setw(8) << std::fixed << std::setprecision(2) << frac.getLength() << " | "
                  << std::setw(8) << std::fixed << std::setprecision(2) << frac.getHeight() << "\n";
    }
    std::cout << std::string(80, '-') << "\n\n";

    double cell_size = computeOptimalCellSize(elements);
    size_t num_fractures = elements.size();

    // Initialize connections list
    std::vector<std::vector<int>> connections(num_fractures);

    // Build spatial grid
    SpatialGrid grid(cell_size);

    for (const auto& element : elements) {
        grid.insert(element);
    }

    // Query and check for intersections
    for (const auto& element_i : elements) {
        auto nearby = grid.queryNearby(element_i);
        for (const auto* element_j : nearby) {
            if (element_i.getID() < element_j->getID()) {
                if (Intersection::fracturesIntersect(element_i, *element_j)) {
                    connections[element_i.getID()].push_back(element_j->getID());
                    connections[element_j->getID()].push_back(element_i.getID());
                }
            }
        }
    }

    // Count total unique intersections
    int total_intersections = 0;
    for (const auto& conn : connections) {
        total_intersections += conn.size();
    }
    total_intersections /= 2; // Each intersection counted twice

    std::cout << "Results:\n";
    std::cout << std::string(60, '-') << "\n";
    std::cout << "Total intersections found: " << total_intersections << "\n\n";

    // Print connection table
    if (total_intersections > 0) {
        std::cout << "Fracture Connection Table:\n";
        std::cout << std::string(40, '-') << "\n";
        std::cout << std::setw(12) << "Fracture" << " | " << "Connects To\n";
        std::cout << std::string(40, '-') << "\n";
        
        for (size_t i = 0; i < num_fractures; ++i) {
            if (!connections[i].empty()) {
                std::cout << std::setw(12) << i << " | ";
                for (size_t j = 0; j < connections[i].size(); ++j) {
                    std::cout << connections[i][j];
                    if (j < connections[i].size() - 1) std::cout << ", ";
                }
                std::cout << "\n";
            }
        }
        std::cout << std::string(40, '-') << "\n\n";
    }

    // Print connections for each fracture
    for (size_t i = 0; i < num_fractures; ++i) {
        elements[i].print();
        if (connections[i].empty()) {
            std::cout << "  No intersections\n";
        } else {
            std::cout << "  Intersects with fracture(s): ";
            for (size_t j = 0; j < connections[i].size(); ++j) {
                std::cout << connections[i][j];
                if (j < connections[i].size() - 1) std::cout << ", ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
    }

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
    
    std::cout << "Cell size: " << cell_size << " m\n";
    std::cout << "Number of occupied cells: " << num_cells << "\n";
    std::cout << "Total fracture references in grid: " << num_fractures_in_grid << "\n";
    std::cout << "Fill ratio (fractures per cell): " << std::fixed << std::setprecision(2) << fill_ratio << "\n";
    std::cout << "Domain bounds:\n";
    std::cout << "  X: [" << std::fixed << std::setprecision(2) << domain_min.x() << ", " << domain_max.x() 
              << "] (width = " << domain_width << " m)\n";
    std::cout << "  Y: [" << domain_min.y() << ", " << domain_max.y() 
              << "] (height = " << domain_height << " m)\n";
    std::cout << "  Z: [" << domain_min.z() << ", " << domain_max.z() 
              << "] (depth = " << domain_depth << " m)\n";
    std::cout << "\n";

    // Draw ASCII map of fractures using actual domain bounds
    double map_width = domain_width * 1.1;
    double map_height = domain_height * 1.1;
    drawASCIIMap(elements, map_width, map_height, 100, 40);
    std::cout << "\n";
    drawSpatialGridMap(grid, cell_size, 50, 40);
}

int main() {
    // Create timestamped output file
    std::time_t now = std::time(nullptr);
    std::tm* local_time = std::localtime(&now);
    char timestamp[32];
    std::strftime(timestamp, sizeof(timestamp), "%Y%m%d_%H%M%S", local_time);
    std::string output_filename = "fracture_output_" + std::string(timestamp) + ".txt";
    
    std::ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        std::cerr << "Warning: Could not create output file '" << output_filename << "'\n";
        std::cerr << "Continuing with console output only.\n\n";
    } else {
        std::cout << "Output will be saved to: " << output_filename << "\n\n";
    }
    
    // Set up tee buffer to write to both console and file
    std::streambuf* cout_backup = std::cout.rdbuf();
    TeeBuf tee_buf(cout_backup, output_file.rdbuf());
    std::cout.rdbuf(&tee_buf);

    // Run multiple examples
    
    // Example 1: Cross-pattern intersections
    runExample(
        "Cross-Pattern Intersections",
        "Four fractures arranged in a cross pattern to demonstrate multiple intersections.",
        {
            -5.0,  0.0,  0.0,   // F0: left side, vertical
            5.0,   0.0,  0.0,   // F1: right side, vertical
            0.0,  -5.0,  0.0,   // F2: bottom, horizontal
            0.0,   5.0,  0.0    // F3: top, horizontal
        },
        {90.0, 90.0, 0.0, 0.0},
        12.0,  // length
        6.0    // height
    );

    // Example 2: Radial fractures
    runExample(
        "Radial Fracture Pattern",
        "Six fractures radiating from a common center point at 60° intervals.",
        {
            0.0,  0.0,  0.0,   // F0: center
            0.0,  0.0,  0.0,   // F1: center
            0.0,  0.0,  0.0,   // F2: center
            0.0,  0.0,  0.0,   // F3: center
            0.0,  0.0,  0.0,   // F4: center
            0.0,  0.0,  0.0    // F5: center
        },
        {0.0, 30.0, 60.0, 90.0, 120.0, 150.0},
        10.0,  // length
        5.0    // height
    );

    // Example 3: Non-intersecting parallel fractures
    runExample(
        "Parallel Non-Intersecting Fractures",
        "Three parallel vertical fractures with no intersections.",
        {
            -10.0,  0.0,  0.0,   // F0: left
            0.0,    0.0,  0.0,   // F1: center
            10.0,   0.0,  0.0    // F2: right
        },
        {90.0, 90.0, 90.0},
        8.0,   // length
        4.0    // height
    );

    // Example 4: Complex 3D intersections
    runExample(
        "3D Multi-Level Fracture Network",
        "Eight fractures distributed across different Z levels forming a complex network.",
        {
            -8.0,  -8.0,  -5.0,  // F0: lower level
            8.0,   -8.0,  -5.0,  // F1: lower level
            -8.0,   8.0,   5.0,  // F2: upper level
            8.0,    8.0,   5.0,  // F3: upper level
            0.0,   -8.0,   0.0,  // F4: mid level
            0.0,    8.0,   0.0,  // F5: mid level
            -8.0,   0.0,   0.0,  // F6: mid level
            8.0,    0.0,   0.0   // F7: mid level
        },
        {45.0, 135.0, 45.0, 135.0, 0.0, 0.0, 90.0, 90.0},
        15.0,  // length
        8.0    // height
    );

    // Example 5: Dense fracture cluster
    runExample(
        "Dense Fracture Cluster",
        "Ten fractures in a small area with varying orientations.",
        {
            -2.0,  0.0,  0.0,
            2.0,   0.0,  0.0,
            0.0,  -2.0,  0.0,
            0.0,   2.0,  0.0,
            -1.5,  1.5,  0.0,
            1.5,   1.5,  0.0,
            -1.5, -1.5,  0.0,
            1.5,  -1.5,  0.0,
            0.0,   0.0,  2.0,
            0.0,   0.0, -2.0
        },
        {0.0, 90.0, 45.0, 135.0, 30.0, 60.0, 120.0, 150.0, 0.0, 90.0},
        6.0,   // length
        3.0    // height
    );

    // Example 6: XY overlap but no Z intersection
    runExample(
        "XY Overlap with Z Separation (No Intersection)",
        "Two fractures with same XY position and orientation but different Z levels - should NOT intersect.",
        {
            0.0,  0.0,  -10.0,   // F0: lower level (z = -10)
            0.0,  0.0,   10.0    // F1: upper level (z = +10)
        },
        {45.0, 45.0},   // Same orientation
        12.0,   // length
        5.0     // height (total Z extent from -12.5 to -7.5 and +7.5 to +12.5)
    );

    // Example 7: Random fractures
    {
        // User-defined parameters for random generation
        int num_fractures = 30;
        double domain_x_min = 0.0, domain_x_max = 500.0;
        double domain_y_min = 0.0, domain_y_max = 500.0;
        double domain_z_min = 0.0, domain_z_max = 500.0;
        double min_length = 5.0, max_length = 50.0;
        double min_height = 5.0, max_height = 50.0; 
        unsigned int seed = 12345;
        
        // Build description with parameters
        std::ostringstream desc;
        desc << "Generate a user-defined number of randomly oriented fractures over a specified domain.\n\n"
             << "Random Generation Parameters:\n"
             << "  Domain: X[" << domain_x_min << ", " << domain_x_max << "]"
             << " Y[" << domain_y_min << ", " << domain_y_max << "]"
             << " Z[" << domain_z_min << ", " << domain_z_max << "]\n"
             << "  Fracture length range: [" << min_length << ", " << max_length << "]\n"
             << "  Fracture height range: [" << min_height << ", " << max_height << "]\n"
             << "  Random seed: " << seed;
        
        // Generate random fractures
        std::vector<FractureElement> fractures = generateRandomFractures(
            num_fractures,
            domain_x_min, domain_x_max,
            domain_y_min, domain_y_max,
            domain_z_min, domain_z_max,
            max_length, max_length,
            max_height, max_height,
            seed
        );
        
        // Run the example using the overloaded function
        runExample("Random Fracture Generation", desc.str(), fractures);
    }
    
    // Restore original cout buffer
    std::cout.rdbuf(cout_backup);
    
    // Close output file
    if (output_file.is_open()) {
        output_file.close();
        std::cout << "\nOutput saved to: " << output_filename << "\n";
    }
    
    return 0;
}
