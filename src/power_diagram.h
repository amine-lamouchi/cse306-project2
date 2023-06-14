#pragma once

#include <vector>
#include "vector.h"
#include "polygon.h"
#include <string>

void save_svg(const std::vector<Polygon> &polygons, std::string filename, std::string fillcol = "none") {
        FILE* f = fopen(filename.c_str(), "w+"); 
        fprintf(f, "<svg xmlns = \"http://www.w3.org/2000/svg\" width = \"1000\" height = \"1000\">\n");
        for (int i=0; i<polygons.size(); i++) {
            fprintf(f, "<g>\n");
            fprintf(f, "<polygon points = \""); 
            for (int j = 0; j < polygons[i].vertices.size(); j++) {
                fprintf(f, "%3.3f, %3.3f ", (polygons[i].vertices[j][0] * 1000), (1000 - polygons[i].vertices[j][1] * 1000));
            }
            fprintf(f, "\"\nfill = \"%s\" stroke = \"black\"/>\n", fillcol.c_str());
            fprintf(f, "</g>\n");
        }
        fprintf(f, "</svg>\n");
        fclose(f);
}
 
class PowerDiagram {
    public:
    PowerDiagram() {
        const int N_disk = 30;
        disk.vertices.resize(N_disk);

        for (int i = 0; i < N_disk; i++) {
            double t = 2 * M_PI * i / N_disk;
            disk.vertices[i][0] = cos(t);
            disk.vertices[i][1] = sin(t);
            disk.vertices[i][2] = 0;
        }
    }

    PowerDiagram(const std::vector<Vector> &pts, const std::vector<double> &wghts) {
        points = pts;
        weights = wghts;
        const int N_disk = 100;
        powerDiagram.resize(points.size());
        disk.vertices.resize(N_disk);

        for (int i = 0; i < N_disk; i++) {
            double t = 2 * M_PI * i / N_disk;
            disk.vertices[i][0] = cos(t);
            disk.vertices[i][1] = sin(t);
            disk.vertices[i][2] = 0;
        }
    }

    Polygon clipPolygonByEdge(const Polygon &poly, const Vector &u, const Vector &v) {
        Polygon result;
        result.vertices.reserve(poly.vertices.size()+1);
        Vector N(v[1]-u[1], u[0]-v[0], 0);

        for (int i = 0; i < poly.vertices.size(); i++) {
            const Vector &A = (i==0) ? poly.vertices[poly.vertices.size()-1]: poly.vertices[i-1];
            const Vector &B = poly.vertices[i];
            double t = dot(u-A, N)/dot(B-A, N);
            Vector P = A + t * (B-A);

            if (dot(B-u, N) < 0) {
                if (dot(A-u, N) > 0) {
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            } else {
                if (dot(A-u, N) < 0) {
                    result.vertices.push_back(P);
                }
            }
        }
        return result;
    }

    Polygon clipPolygonByBissector(const Polygon &poly, int index_0, int index_i, const Vector &P0, const Vector &Pi) {
        Vector M = (P0 + Pi) / 2.0;
        Vector Mprime = M + (weights[index_0] - weights[index_i])/(2.*(P0-Pi).norm2()) * (Pi-P0);
        Polygon result;

        for (int i = 0; i < poly.vertices.size(); i++) {
            const Vector &A = (i==0) ? poly.vertices[poly.vertices.size()-1]: poly.vertices[i-1];
            const Vector &B = poly.vertices[i];
            double t = dot(Mprime-A, Pi-P0)/dot(B-A, Pi - P0);
            Vector P = A + t * (B-A);

            if ((B-P0).norm2() - weights[index_0] < (B-Pi).norm2() - weights[index_i]) {
                if ((A-P0).norm2() - weights[index_0] > (A-Pi).norm2() - weights[index_i]) {
                    result.vertices.push_back(P);
                }
                result.vertices.push_back(B);
            } else {
                if ((A-P0).norm2() - weights[index_0] < (A-Pi).norm2() - weights[index_i]) {
                    result.vertices.push_back(P);
                }
            }
        }
        return result;
    }

    Polygon intersect_with_disk(const Polygon& polygon, const Vector& center, double radius) {
        Polygon result(polygon);
        for (int i=0; i<disk.vertices.size(); i++) {
            const Vector& u = disk.vertices[i]*radius+center;
            const Vector& v = disk.vertices[(i+1)%30]*radius+center;
            result = clipPolygonByEdge(result, u, v);
        }
        return result;
    }

    Polygon computePowerCell(int idx) {
        Polygon result;
        result.vertices.resize(4);
        result.vertices[0] = Vector(0, 0, 0);
        result.vertices[1] = Vector(0, 1, 0);
        result.vertices[2] = Vector(1, 1, 0);
        result.vertices[3] = Vector(1, 0, 0);

        for (int i = 0; i < points.size(); i++) {
            if (i == idx) continue;
            result = clipPolygonByBissector(result, idx, i, points[idx], points[i]);
        }
        result = intersect_with_disk(result, points[idx], std::sqrt(weights[idx]-weights[weights.size()-1]));
        return result;
    }

    void compute() {
        powerDiagram.resize(points.size());

        for (int i = 0; i < points.size(); i++) {
            powerDiagram[i] = computePowerCell(i);
        }
    }

    void save(std::string filename) const {
        save_svg(powerDiagram, filename, "blue");
    }


    std::vector<Polygon> powerDiagram;
    std::vector<Vector> points;
    std::vector<double> weights;
    Polygon disk;
};
