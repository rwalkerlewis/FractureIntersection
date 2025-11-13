#ifndef GEOMETRY_H
#define GEOMETRY_H

#include <cmath>
#include <array>

// ========================================================================= 80
// Generic Point Class (2D/3D)
// ========================================================================= 80
/*
* Generic point class that works for 2D and 3D points.
* Template parameter Dim specifies dimensionality (2 or 3).
*/
template<int Dim>
class Point {
private:
    std::array<double, Dim> coords_;

public:
    // Default constructor
    Point() {
        coords_.fill(0.0);
    }

    // Constructor for 2D points
    template<int D = Dim, typename = std::enable_if_t<D == 2>>
    Point(double x_, double y_) {
        coords_[0] = x_;
        coords_[1] = y_;
    }

    // Constructor for 3D points
    template<int D = Dim, typename = std::enable_if_t<D == 3>>
    Point(double x_, double y_, double z_) {
        coords_[0] = x_;
        coords_[1] = y_;
        coords_[2] = z_;
    }

    // Property accessors for backward compatibility
    double& x() { return coords_[0]; }
    const double& x() const { return coords_[0]; }
    
    double& y() { return coords_[1]; }
    const double& y() const { return coords_[1]; }
    
    template<int D = Dim, typename = std::enable_if_t<D == 3>>
    double& z() { return coords_[2]; }
    
    template<int D = Dim, typename = std::enable_if_t<D == 3>>
    const double& z() const { return coords_[2]; }

    // Internal coords array access
    std::array<double, Dim>& coords() { return coords_; }
    const std::array<double, Dim>& coords() const { return coords_; }

    // Operator overloads
    Point operator-(const Point& other) const {
        Point result;
        for (int i = 0; i < Dim; ++i) {
            result.coords_[i] = coords_[i] - other.coords_[i];
        }
        return result;
    }

    Point operator+(const Point& other) const {
        Point result;
        for (int i = 0; i < Dim; ++i) {
            result.coords_[i] = coords_[i] + other.coords_[i];
        }
        return result;
    }

    Point operator*(double scalar) const {
        Point result;
        for (int i = 0; i < Dim; ++i) {
            result.coords_[i] = coords_[i] * scalar;
        }
        return result;
    }

    double dot(const Point& other) const {
        double result = 0.0;
        for (int i = 0; i < Dim; ++i) {
            result += coords_[i] * other.coords_[i];
        }
        return result;
    }

    double length() const {
        return std::sqrt(dot(*this));
    }

    // Cross product (3D only)
    template<int D = Dim, typename = std::enable_if_t<D == 3>>
    Point cross(const Point& other) const {
        Point result;
        result.coords_[0] = coords_[1] * other.coords_[2] - coords_[2] * other.coords_[1];
        result.coords_[1] = coords_[2] * other.coords_[0] - coords_[0] * other.coords_[2];
        result.coords_[2] = coords_[0] * other.coords_[1] - coords_[1] * other.coords_[0];
        return result;
    }

    // 2D cross product (returns scalar)
    template<int D = Dim, typename = std::enable_if_t<D == 2>>
    double cross(const Point& other) const {
        return coords_[0] * other.coords_[1] - coords_[1] * other.coords_[0];
    }
};

/*
* Axis-aligned bounding box for a fracture element
*/
struct BoundingBox {
    Point<3> min, max;

    /*
    * Check if this bounding box intersects with another bounding box
    */
    /* 
    * Check to see if this bounding box intersects with another bounding box
    */
    bool intersects(const BoundingBox& other) const {
        return (max.x() >= other.min.x() && min.x() <= other.max.x()) &&
               (max.y() >= other.min.y() && min.y() <= other.max.y()) &&
               (max.z() >= other.min.z() && min.z() <= other.max.z());
    }

    /*
    * Check if this bounding box intersects with another bounding box
    */
    void expand(const Point<3>& point) {
        min.x() = std::min(min.x(), point.x());
        min.y() = std::min(min.y(), point.y());
        min.z() = std::min(min.z(), point.z());
        max.x() = std::max(max.x(), point.x());
        max.y() = std::max(max.y(), point.y());
        max.z() = std::max(max.z(), point.z());
    }

    /* 
    * Get center of bounding box
    */
    Point<3> center() const {
        return Point<3>(
            (min.x() + max.x()) / 2.0,
            (min.y() + max.y()) / 2.0,
            (min.z() + max.z()) / 2.0);
    }

};

#endif // GEOMETRY_H