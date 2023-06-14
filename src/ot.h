#pragma once

#include <vector>
#include "vector.h"
#include "power_diagram.h"
#include "lbfgs.h"
#include <iostream>
#define VOLUME_AIR 0.7
#define VOLUME_FLUID 0.3

class OT {
public:
    OT() {};    
    OT(std::vector<Vector> &points, const std::vector<double> &lambdas) {
        this->points = points;
        this->lambdas = lambdas;
    };

    void solve() {
        solution.points = points;
        solution.weights.resize(points.size()+1);
        std::fill(solution.weights.begin(), solution.weights.end(), 1.0);
        solution.weights[solution.weights.size()-1] = 1- 0.001;

        double fx = 0;

        int ret = lbfgs(points.size()+1, &solution.weights[0], &fx, _evaluate, _progress, this, NULL);
        solution.compute();
    }

    static lbfgsfloatval_t _evaluate(
        void *instance,
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        return reinterpret_cast<OT*>(instance)->evaluate(x, g, n, step);
    }

    lbfgsfloatval_t evaluate(
        const lbfgsfloatval_t *x,
        lbfgsfloatval_t *g,
        const int n,
        const lbfgsfloatval_t step
        )
    {
        lbfgsfloatval_t fx = 0.0;
    
        for (int i = 0; i < n; i++)
            solution.weights[i] = x[i];
        solution.compute();

        double s1 = 0;
        double s2 = 0;
        double s3 = 0;
        double estimated_volume_fluid = 0;
        for(int i = 0; i < n-1; ++i) {
            g[i] = -(lambdas[i] - solution.powerDiagram[i].area());
            s3 += lambdas[i] * x[i];
            s2 -= x[i] * solution.powerDiagram[i].area();
            s1 += solution.powerDiagram[i].integrateSquareDistance(solution.points[i]);
            estimated_volume_fluid += solution.powerDiagram[i].area();
        }

        double estimated_volume_air = 1.0 - estimated_volume_fluid;
        fx = s1 + s2 + s3 + x[n-1] * (VOLUME_AIR - estimated_volume_air);
        g[n-1] = estimated_volume_air - VOLUME_AIR;
        return -fx;
    }

    static int _progress(
        void *instance,
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        return reinterpret_cast<OT*>(instance)->progress(x, g, fx, xnorm, gnorm, step, n, k, ls);
    }

    int progress(
        const lbfgsfloatval_t *x,
        const lbfgsfloatval_t *g,
        const lbfgsfloatval_t fx,
        const lbfgsfloatval_t xnorm,
        const lbfgsfloatval_t gnorm,
        const lbfgsfloatval_t step,
        int n,
        int k,
        int ls
        )
    {
        for (int i = 0; i < n; i++) {
            solution.weights[i] = x[i];
        }
        solution.compute();

        double max_diff = 0;
        for (int i = 0; i < n-1; i++) {
            max_diff = std::max(max_diff, std::abs(solution.powerDiagram[i].area() - lambdas[i]));
        }

        std::cout << "fx: " << fx << " max_diff= " << max_diff << "\t gnorm= " << gnorm << std::endl;



        printf("Iteration %d:\n", k);
        printf("  fx = %f, x[0] = %f, x[1] = %f\n", fx, x[0], x[1]);
        printf("  xnorm = %f, gnorm = %f, step = %f\n", xnorm, gnorm, step);
        printf("\n");
        return 0;
    }

    std::vector<Vector> points;
    std::vector<double> lambdas;
    PowerDiagram solution;
};