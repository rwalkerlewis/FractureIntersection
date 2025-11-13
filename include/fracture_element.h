#ifndef FRACTURE_ELEMENT_H
#define FRACTURE_ELEMENT_H

#include "geometry.h"
#include <cmath>
#include <iostream>
#include <iomanip>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

// ========================================================================= 80
// Fracture Element
// ========================================================================= 80
class FractureElement {
public: 
    // Constructor without dip parameter
    FractureElement(const Point<3>& center, double theta, double length, double height, int id)
        : center_(center), theta_(theta), length_(length), height_(height), id_(id) {
        computeBoundingBox();
    }

    // Data Retrieval
    Point<3> getCenter() const {return center_;}
    double getTheta() const {return theta_;}
    double getLength() const {return length_;}
    double getHeight() const {return height_;}
    int getID() const {return id_;}
    const BoundingBox& getBoundingBox() const {return bounding_box_;}

    /**
     * Return vertical range of the fracture element as a pair of (min_z, max_z)
     */
    void getVerticalRange(double& z_min, double& z_max) const {
        double half_height = height_ / 2.0;
        z_min = center_.z() - half_height;
        z_max = center_.z() + half_height;
    }

    /**
     * Get the two endpoints of the fracture's horizontal projection
     * The fracture projects to a line segment in the (x,y) plane
     */
    void getHorizontalEndpoints(Point<2>& p1, Point<2>& p2) const {
        double half_length = length_ / 2.0;
        double theta_rad = theta_ * M_PI / 180.0; // Convert to radians

        double dx = half_length * cos(theta_rad);
        double dy = half_length * sin(theta_rad);

        p1.x() = center_.x() - dx;
        p1.y() = center_.y() - dy;
        p2.x() = center_.x() + dx;
        p2.y() = center_.y() + dy;
    }

    /*
    * Compute Bounding Box for this fracture element
    */
    void computeBoundingBox() {
        Point<2> p1, p2;
        getHorizontalEndpoints(p1, p2);
        
        double z_min, z_max;
        getVerticalRange(z_min, z_max);

        bounding_box_.min.x() = std::min(p1.x(), p2.x());
        bounding_box_.min.y() = std::min(p1.y(), p2.y());
        bounding_box_.min.z() = z_min;

        bounding_box_.max.x() = std::max(p1.x(), p2.x());
        bounding_box_.max.y() = std::max(p1.y(), p2.y());
        bounding_box_.max.z() = z_max;
    }

    // Print methods
    void print() const {
        std::cout << "FractureElement " << id_ << ":\n"
                  << "  Center: (" << std::fixed << std::setprecision(2) 
                  << center_.x() << ", " << center_.y() << ", " << center_.z() << ")\n"
                  << "  Theta (orientation): " << theta_ << " degrees\n"
                  << "  Length: " << length_ << " m\n"
                  << "  Height: " << height_ << " m\n";
    }
    
    friend std::ostream& operator<<(std::ostream& os, const FractureElement& element);

private:
    Point<3> center_;
    double theta_;
    double length_;
    double height_;
    int id_;
    BoundingBox bounding_box_;
};

// Overload << operator for easy printing
inline std::ostream& operator<<(std::ostream& os, const FractureElement& element) {
    os << "FractureElement " << element.getID() << " ["
       << "center=(" << std::fixed << std::setprecision(2)
       << element.getCenter().x() << ", "
       << element.getCenter().y() << ", "
       << element.getCenter().z() << "), "
       << "theta=" << element.getTheta() << "°, "
       << "length=" << element.getLength() << "m, "
       << "height=" << element.getHeight() << "m]";
    return os;
}

#endif // FRACTURE_ELEMENT_H