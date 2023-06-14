#include <vector>
#include "vector.h"
#include "fluid.h"
#include "ot.h"
#include <cassert>
#include "polygon.h"
#include <chrono>

int main() {
    Fluid fluid(50);
    auto start = std::chrono::high_resolution_clock::now();
    fluid.runFluid();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
    std::cout << "Time: " << duration.count() << " seconds" << std::endl;
    std::vector<Vector> points(32);
    std::vector<double> lambdas(32);

    for (int i = 0; i < points.size(); i++) {
        points[i][0] = rand() / static_cast<double>(RAND_MAX);
        points[i][1] = rand() / static_cast<double>(RAND_MAX);
        points[i][2] = 0;
        lambdas[i] = 1.0 / points.size();
    }

    OT ot(points, lambdas);
    ot.solve();
    ot.solution.save("voronoi.svg");
    return 0;
}