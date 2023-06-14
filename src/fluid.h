#pragma once

#include <vector>
#include "vector.h"
#include <string>
#include "polygon.h"
#include <omp.h>
#include <sstream>
#include "stb_image_write.h"
#include "ot.h"

void save_frame(const std::vector<Polygon>& cells, std::string filename, int N, int frameid = 0) {
	int W = 512, H = 512;
	std::vector<unsigned char> image(W * H * 3, 255);
#pragma omp parallel for schedule(dynamic)
	for (int i = 0; i < cells.size(); i++) {

		double bminx = 1E9, bminy = 1E9, bmaxx = -1E9, bmaxy = -1E9;
		for (int j = 0; j < cells[i].vertices.size(); j++) {
			bminx = std::min(bminx, cells[i].vertices[j][0]);
			bminy = std::min(bminy, cells[i].vertices[j][1]);
			bmaxx = std::max(bmaxx, cells[i].vertices[j][0]);
			bmaxy = std::max(bmaxy, cells[i].vertices[j][1]);
		}
		bminx = std::min(W - 1., std::max(0., W * bminx));
		bminy = std::min(H - 1., std::max(0., H * bminy));
		bmaxx = std::max(W - 1., std::max(0., W * bmaxx));
		bmaxy = std::max(H - 1., std::max(0., H * bmaxy));

		for (int y = bminy; y < bmaxy; y++) {
			for (int x = bminx; x < bmaxx; x++) {
				int prevSign = 0;
				bool isInside = true;
				double mindistEdge = 1E9;
				for (int j = 0; j < cells[i].vertices.size(); j++) {
					double x0 = cells[i].vertices[j][0] * W;
					double y0 = cells[i].vertices[j][1] * H;
					double x1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][0] * W;
					double y1 = cells[i].vertices[(j + 1) % cells[i].vertices.size()][1] * H;
					double det = (x - x0) * (y1 - y0) - (y - y0) * (x1 - x0);
					int sign = (det > 0) ? 1 : ((det < 0) ? -1 : 0);
					if (prevSign == 0) prevSign = sign; else
						if (sign == 0) sign = prevSign; else
							if (sign != prevSign) {
								isInside = false;
								break;
							}
					prevSign = sign;
					double edgeLen = sqrt((x1 - x0) * (x1 - x0) + (y1 - y0) * (y1 - y0));
					double distEdge = std::abs(det) / edgeLen;
					double dotp = (x - x0) * (x1 - x0) + (y - y0) * (y1 - y0);
					if (dotp<0 || dotp>edgeLen * edgeLen) distEdge = 1E9;
					mindistEdge = std::min(mindistEdge, distEdge);
				}
				if (isInside) {
					{   // the N first particles may represent fluid, displayed in blue
						image[((H - y - 1)*W + x) * 3] = 0;
						image[((H - y - 1)*W + x) * 3 + 1] = 0;
						image[((H - y - 1)*W + x) * 3 + 2] = 255;
					}
					if (mindistEdge <= 2) {
						image[((H - y - 1) * W + x) * 3] = 0;
						image[((H - y - 1) * W + x) * 3 + 1] = 0;
						image[((H - y - 1) * W + x) * 3 + 2] = 0;
					}

				}

			}
		}
	}
	std::ostringstream os;
	os << filename << frameid << ".png";
	stbi_write_png(os.str().c_str(), W, H, 3, &image[0], 0);
}

class Fluid
{
public:

    Fluid(int N)
    {
        this->N = N;
        particles.resize(N);
        for (int i = 0; i < N; i++)
        {
            particles[i] = Vector(rand() / (double)RAND_MAX, rand() / (double)RAND_MAX, 0.);
        }
        velocities.resize(N, Vector(0., 0., 0.));
    };

    void stepFluid()
    {
        otsolver.points = particles;
        otsolver.lambdas = std::vector<double>(particles.size(), 1.0 / particles.size() * VOLUME_FLUID);
        otsolver.solve();
        const double dt = 0.002;
        double mass_particle = 50;
        double epsilon2 = 0.004 * 0.004;
        for (int i = 0; i < particles.size(); i++)
        {
            Vector gravity = Vector(0, -9.81, 0) * mass_particle;
            Vector centroid = otsolver.solution.powerDiagram[i].centroid();
            Vector otForce = 1. / epsilon2 * (centroid - particles[i]);
            Vector forces = gravity + otForce;
            velocities[i] += dt / mass_particle * forces;
            particles[i] += dt * velocities[i];
        }
    }

    void runFluid()
    {
        for (int i = 0; i < N; i++)
        {
            stepFluid();
            save_frame(otsolver.solution.powerDiagram, "./animation_frames/animation", N, i);
        }
    }
    OT otsolver;
    std::vector<Vector> particles;
    std::vector<Vector> velocities;
    int N;
};