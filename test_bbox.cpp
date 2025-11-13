#include <iostream>
#include <cmath>
#include "fracture_element.h"

int main() {
    Point<3> center(-5.0, 0.0, 0.0);
    FractureElement elem(center, 90.0, 12.0, 6.0, 0);
    
    auto bbox = elem.getBoundingBox();
    std::cout << "BBox min: (" << bbox.min.x() << ", " << bbox.min.y() << ", " << bbox.min.z() << ")\n";
    std::cout << "BBox max: (" << bbox.max.x() << ", " << bbox.max.y() << ", " << bbox.max.z() << ")\n";
    
    return 0;
}
