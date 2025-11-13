#ifndef SPATIAL_GRID_H
#define SPATIAL_GRID_H 

#include "fracture_element.h"
#include <vector>
#include <cmath>
#include <unordered_map>
#include <unordered_set>

/**
 * Spatial Grid for efficient fracture intersection queries
 * 
 * KEY EFFICIENCY STRATEGY:
 * - Uses hash map (unordered_map) instead of dense array
 * - Only stores grid cells that CONTAIN fractures (sparse representation)
 * - Automatically adapts to domain size without wasting empty cells
 * - Tracks world-space bounds to limit grid to fracture domain
 * - Fill ratio < 1.0 indicates sparse usage (good efficiency)
 */
class SpatialGrid {
private:
    double cell_size_;
    // Bounds of the domain containing fractures (in world coordinates)
    Point<3> domain_min_, domain_max_;
    bool bounds_initialized_;

public:
    explicit SpatialGrid(double cell_size) : cell_size_(cell_size), bounds_initialized_(false) {
        if (cell_size <= 0.0) {
            throw std::invalid_argument("Cell size must be positive");
        }
        domain_min_ = Point<3>(0.0, 0.0, 0.0);
        domain_max_ = Point<3>(0.0, 0.0, 0.0);
    }

    void insert(const FractureElement& element) {
        const BoundingBox& bbox = element.getBoundingBox();
        
        // Update domain bounds
        if (!bounds_initialized_) {
            domain_min_ = bbox.min;
            domain_max_ = bbox.max;
            bounds_initialized_ = true;
        } else {
            // Expand bounds to include this element
            domain_min_.x() = std::min(domain_min_.x(), bbox.min.x());
            domain_min_.y() = std::min(domain_min_.y(), bbox.min.y());
            domain_min_.z() = std::min(domain_min_.z(), bbox.min.z());
            domain_max_.x() = std::max(domain_max_.x(), bbox.max.x());
            domain_max_.y() = std::max(domain_max_.y(), bbox.max.y());
            domain_max_.z() = std::max(domain_max_.z(), bbox.max.z());
        }
        
        std::vector<GridCell> cells = getGridCellsForBoundingBox(bbox);

        for (const auto& cell : cells) {
            grid_[cell].push_back(&element);
        }
    }

    std::vector<const FractureElement*> queryNearby(const FractureElement& fracture) const {
        // Find all grid cell fracture overlaps
        std::vector<GridCell> cells = getGridCellsForBoundingBox(fracture.getBoundingBox());

        // Use a set to avoid duplicates
        std::unordered_set<const FractureElement*> nearby_set;

        for (const auto& cell : cells) {
            auto it = grid_.find(cell);
            if (it != grid_.end()) {
                for (const auto* element : it->second) {
                    // Do not include the fracture itself
                    if (element->getID() != fracture.getID()) {
                        nearby_set.insert(element);
                    }
                }
            }
        }
        return std::vector<const FractureElement*>(nearby_set.begin(), nearby_set.end());
    }
    void clear() {
        grid_.clear();
    }

    double getCellSize() const {
        return cell_size_;
    }

    // Get the grid bounds (min and max cell coordinates)
    void getGridBounds(int& min_x, int& max_x, int& min_y, int& max_y, int& min_z, int& max_z) const {
        min_x = max_x = min_y = max_y = min_z = max_z = 0;
        bool first = true;
        for (const auto& pair : grid_) {
            const GridCell& cell = pair.first;
            if (first) {
                min_x = max_x = cell.x;
                min_y = max_y = cell.y;
                min_z = max_z = cell.z;
                first = false;
            } else {
                min_x = std::min(min_x, cell.x);
                max_x = std::max(max_x, cell.x);
                min_y = std::min(min_y, cell.y);
                max_y = std::max(max_y, cell.y);
                min_z = std::min(min_z, cell.z);
                max_z = std::max(max_z, cell.z);
            }
        }
    }

    // Get fractures in a specific cell
    std::vector<int> getFracturesInCell(int x, int y, int z) const {
        std::vector<int> result;
        GridCell cell(x, y, z);
        auto it = grid_.find(cell);
        if (it != grid_.end()) {
            for (const auto* elem : it->second) {
                result.push_back(elem->getID());
            }
        }
        return result;
    }

    // Get the world-space domain bounds (min/max of all fractures)
    void getDomainBounds(Point<3>& min_pt, Point<3>& max_pt) const {
        min_pt = domain_min_;
        max_pt = domain_max_;
    }

    // Get statistics about the grid
    void getGridStatistics(int& num_cells, int& num_fractures, double& fill_ratio) const {
        num_cells = grid_.size();
        std::unordered_set<int> unique_fractures;
        for (const auto& pair : grid_) {
            for (const auto* elem : pair.second) {
                unique_fractures.insert(elem->getID());
            }
        }
        num_fractures = unique_fractures.size();
        
        // Compute theoretical grid size (grid cells needed to cover domain)
        if (bounds_initialized_ && cell_size_ > 0.0) {
            int x_cells = static_cast<int>(std::ceil((domain_max_.x() - domain_min_.x()) / cell_size_)) + 1;
            int y_cells = static_cast<int>(std::ceil((domain_max_.y() - domain_min_.y()) / cell_size_)) + 1;
            int z_cells = static_cast<int>(std::ceil((domain_max_.z() - domain_min_.z()) / cell_size_)) + 1;
            int theoretical_cells = std::max(1, x_cells * y_cells * z_cells);
            fill_ratio = static_cast<double>(num_cells) / theoretical_cells;
        } else {
            fill_ratio = 0.0;
        }
    }

private:

    struct GridCell {
        int x, y, z;

        GridCell(int x_, int y_, int z_) : x(x_), y(y_), z(z_) {}

        bool operator==(const GridCell& other) const {
            return x == other.x && y == other.y && z == other.z;
        }
    };

    struct GridCellHash {
        std::size_t operator()(const GridCell& cell) const {
            std::size_t h1 = std::hash<int>()(cell.x);
            std::size_t h2 = std::hash<int>()(cell.y);
            std::size_t h3 = std::hash<int>()(cell.z);
            return h1 ^ (h2 << 1) ^ (h3 << 2);
        }
    };

    std::unordered_map<GridCell, std::vector<const FractureElement*>, GridCellHash> grid_;

    GridCell pointToCell(const Point<3>& point) const {
        return GridCell(static_cast<int>(std::floor(point.x() / cell_size_)),
                        static_cast<int>(std::floor(point.y() / cell_size_)),
                        static_cast<int>(std::floor(point.z() / cell_size_)));
    }

    std::vector<GridCell> getGridCellsForBoundingBox(const BoundingBox& bbox) const {
        std::vector<GridCell> cells;

        int x_min = static_cast<int>(std::floor(bbox.min.x() / cell_size_));
        int x_max = static_cast<int>(std::floor(bbox.max.x() / cell_size_));
        int y_min = static_cast<int>(std::floor(bbox.min.y() / cell_size_));
        int y_max = static_cast<int>(std::floor(bbox.max.y() / cell_size_));
        int z_min = static_cast<int>(std::floor(bbox.min.z() / cell_size_));
        int z_max = static_cast<int>(std::floor(bbox.max.z() / cell_size_));

        for (int x = x_min; x <= x_max; ++x) {
            for (int y = y_min; y <= y_max; ++y) {
                for (int z = z_min; z <= z_max; ++z) {
                    cells.emplace_back(x, y, z);
                }
            }
        }

        return cells;
    }
};

#endif // SPATIAL_GRID_H
