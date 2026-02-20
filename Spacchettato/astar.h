#ifndef ASTAR_H
#define ASTAR_H

#include <stdint.h>
#include <stdbool.h>

#include "graph.h"

typedef struct {
    double dist_cost;
    double time_cost;
    double bound_path;
    uint32_t *path;
    uint32_t path_len;
} AStarResult;

void free_astar_result(AStarResult *r);
AStarResult copy_res(const AStarResult *src);
AStarResult a_star_with_bound(Graph *g, uint32_t start,
                              uint32_t goal, double W, double dist_best, bool foward,
                              uint32_t s, uint32_t d, double budget, double best_feas_sol,
                              double lambda, uint32_t destination);

#endif
