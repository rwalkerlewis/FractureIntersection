#ifndef INTERSECTION_H
#define INTERSECTION_H

#include "fracture_element.h"
#include "geometry.h"

namespace Intersection {

    /**
     * Fast test to check if bounding boxes of two fracture elements intersect.
     */

     /**
      * Test to see if two line segments in 2D intersect. This is used to check if the horizontal projections of two fracture elements intersect.
      */
     bool segmentsIntersect(const Point<2>& p1, const Point<2>& p2, const Point<2>& q1, const Point<2>& q2) {
        // Compute the direction vectors
        Point<2> r = p2 - p1;
        Point<2> s = q2 - q1;

        // Compute the cross product of the direction vectors
        double r_cross_s = r.cross(s);
        Point<2> qp = q1 - p1;
        double cross_qp_r = qp.cross(r);

        const double epsilon = 1e-8; // Tolerance for floating-point comparisons

        // Check if parallel or colinear
        if (std::abs(r_cross_s) < epsilon) {
            if (std::abs(cross_qp_r) < epsilon) {
                // Colinear case: Check if the segments overlap
                double t0 = qp.dot(r) / r.dot(r);
                double t1 = t0 + s.dot(r) / r.dot(r);
         
                if (t0 > t1) {
                    std::swap(t0, t1);
                }

                return (t1 >= 0 && t0 <= 1);
            }
            else {
                return false; // Parallel but not colinear
            }
        }

        // segments not parallel
        double t = qp.cross(s) / r_cross_s;
        double u = qp.cross(r) / r_cross_s; 

        return (t >= 0 && t <= 1) && (u >= 0 && u <= 1);
    }

    /*
     * Test if two fracture elements intersect by checking first vertical range and then horizontal projection.
     */
    bool fracturesIntersect(const FractureElement& element1, const FractureElement& element2) {
        // Check vertical range overlap
        double z_min1, z_max1, z_min2, z_max2;
        element1.getVerticalRange(z_min1, z_max1);
        element2.getVerticalRange(z_min2, z_max2);

        if (z_max1 < z_min2 || z_max2 < z_min1) {
            return false; // No vertical overlap
        }

        // Check horizontal projection overlap
        Point<2> p1_1, p1_2, p2_1, p2_2;
        element1.getHorizontalEndpoints(p1_1, p1_2);
        element2.getHorizontalEndpoints(p2_1, p2_2);

        return segmentsIntersect(p1_1, p1_2, p2_1, p2_2);
    }

}

#endif // INTERSECTION_H