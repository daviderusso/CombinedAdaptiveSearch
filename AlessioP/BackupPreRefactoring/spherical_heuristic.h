//
// Created by alessio on 25/05/25.
//

#ifndef SPHERICAL_HEURISTIC_H
#define SPHERICAL_HEURISTIC_H

#include <stdint.h>
#include <math.h>

//NOTE: The input coordinates are in nanodegrees (to transform them in degrees they must be multiplied by 1,000,000)
// Length data are given in decimeters (0.1m)

#define Pi 0.000000017453292519943295 // Pi / 180 / 1,000,000
#define COEF 61161607.40544 // = Earth radius M * 0.96 * 10


#ifndef DEG_TO_RAD
#define DEG_TO_RAD(deg) ((deg) * (M_PI / 180.0))
#endif

typedef double cost_t;

// Evaluate the distance between two points from their coordinates
cost_t spherical_heuristic_get_explicit(Graph* g, int u, int dest) {

    int32_t lat1 = g->nodes[u].y;
    int32_t lon1 = g->nodes[u].x;
    int32_t lat2 = g->nodes[dest].y;
    int32_t lon2 = g->nodes[dest].x;

    double Delta_lat_squared = pow(fabs(lat1 - lat2) * Pi, 2);
    double cos_mean_lat = cos((lat1 + lat2) * Pi / 2);
    double x = cos_mean_lat * fabs(lon1 - lon2) * Pi;
    return floor(COEF * sqrt(Delta_lat_squared + pow(x, 2)));
}

#endif //SPHERICAL_HEURISTIC_H
