#pragma once

#include <vector>
#include "vector.h"

class Polygon {  
public:

    double area() const {
        double A = 0;
        int N = vertices.size();
        if (N < 3) return 0;
        for (int i = 0; i < N; i++) {
            A += vertices[i][0] * vertices[(i+1)%N][1] - vertices[(i+1)%N][0] * vertices[i][1];
        }
        return std::abs(A/2.);
    }

    double integrateSquareDistance(const Vector& P) const {
        double value = 0;
        if (vertices.size() < 3) return 0;
        for (int i = 1; i < vertices.size() - 1; i++) {
            Vector triangle[3] {vertices[0], vertices[i], vertices[i+1]};

            double local_value = 0;
            for (int k = 0; k < 3; k++) {
                for (int l = k; l < 3; l++) {
                    local_value += dot(triangle[k]-P, triangle[l]-P);
                }
            }
            Vector e1 = triangle[1] - triangle[0];
            Vector e2 = triangle[2] - triangle[0];
            double area_triangle = 0.5*std::abs(e1[1]*e2[0] - e1[0]*e2[1]);
            value += local_value / 6.*area_triangle;
        }
        return value;

    }

    Vector centroid() const {
        Vector c = Vector(0, 0, 0);
        int N = vertices.size();
        double a = area();
        for (int i = 0; i < N; i++) {
            c[0] += -(vertices[i][0] +  vertices[(i+1)%N][0])* (vertices[i][0] * vertices[(i+1)%N][1] - vertices[(i+1)%N][0] * vertices[i][1]);
            c[1] += (vertices[i][1] +  vertices[(i+1)%N][1])* (vertices[i][0] * vertices[(i+1)%N][1] - vertices[(i+1)%N][0] * vertices[i][1]);
        }
        return c/(6*a);
    }
    std::vector<Vector> vertices;
}; 