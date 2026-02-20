#include "astar.h"

#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <stdio.h>

#include "heap.h"

void free_astar_result(AStarResult *r) {
    if (r->path){
        free(r->path);
        r->path = NULL;
    }
}

AStarResult copy_res(const AStarResult *src) {
    AStarResult dst;

    dst.path_len = src->path_len;
    dst.dist_cost = src->dist_cost;
    dst.time_cost = src->time_cost;

    if (src->path != NULL && src->path_len > 0) {
        dst.path = malloc(src->path_len * sizeof(uint32_t));
        if (!dst.path) {
            perror("malloc failed in copy_res");
            exit(EXIT_FAILURE);
        }
        memcpy(dst.path, src->path, src->path_len * sizeof(uint32_t));
    } else {
        dst.path = NULL;
    }

    return dst;
}

AStarResult a_star_with_bound(Graph *g, uint32_t start,
                              uint32_t goal, double W, double dist_best, bool foward,
                              uint32_t s, uint32_t d, double budget, double best_feas_sol,
                              double lambda, uint32_t destination) {
    uint32_t n = g->n_nodes;
    double *cost = malloc(sizeof(double) * n);
    uint32_t *parent = malloc(sizeof(uint32_t) * n);
    uint32_t remaining = 0;
    for (uint32_t i = 0; i < n; i++) {
        cost[i] = DBL_MAX;
        parent[i] = UINT32_MAX;
        (void)remaining;
    }

    MinHeap *heap = createHeap(n);
    cost[start] = 0;
    push(heap, start, 0);

    while (!isEmpty(heap)) {
        HeapNode curr = pop(heap);
        uint32_t u = curr.node;

        if (!g->nodes[u].active) continue;

        if (u == goal) break;

        double current_cost = cost[u];
        double heuristic = 0;

        if (lambda == 0 && foward) g->nodes[u].time_from_source = cost[u];
        if (lambda == 0 && !foward) g->nodes[u].time_to_dest = cost[u];
        if (lambda == 1 && foward) g->nodes[u].dist_from_source = cost[u];
        if (lambda == 1 && !foward) g->nodes[u].dist_to_dest = cost[u];

        if ((g->nodes[u].dist_from_source + g->nodes[u].dist_to_dest >= dist_best) ||
            (g->nodes[u].time_from_source + g->nodes[u].time_to_dest > W)) {
            g->nodes[u].active = false;
            continue;
        }

        if (foward) heuristic = lambda * g->nodes[u].dist_to_dest + (1 - lambda) * g->nodes[u].time_to_dest;
        else heuristic = lambda * g->nodes[u].dist_from_source + (1 - lambda) * g->nodes[u].time_from_source;

        double estimated_total = current_cost + heuristic;

        if (lambda == 0 && estimated_total > W){
            g->nodes[u].active = false;
            continue;
        }

        if (lambda == 1 && estimated_total >= dist_best){
            g->nodes[u].active = false;
            continue;
        }

        if (lambda > 0 && lambda < 1 && estimated_total > (lambda * dist_best + (1- lambda) * W))
            continue;

        for (uint32_t i = 0; i < (foward ? g->nodes[u].out_count : g->nodes[u].in_count); i++) {
            uint32_t v;

            double edge_cost;
            if (foward) {
                Edge e = g->nodes[u].out_edges[i];
                v = e.to;
                edge_cost = lambda * e.distance + (1 - lambda) * e.time;
            } else {
                InEdge e = g->nodes[u].in_edges[i];
                v = e.from;
                edge_cost = lambda * e.distance + (1 - lambda) * e.time;
            }

            if (!g->nodes[v].active) continue;

            double new_cost = cost[u] + edge_cost;

            if (new_cost < cost[v]) {
                cost[v] = new_cost;
                parent[v] = u;

                if (lambda == 0 && foward) g->nodes[v].time_from_source = cost[v];
                if (lambda == 0 && !foward) g->nodes[v].time_to_dest = cost[v];
                if (lambda == 1 && foward) g->nodes[v].dist_from_source = cost[v];
                if (lambda == 1 && !foward) g->nodes[v].dist_to_dest = cost[v];

                double h = 0;

                if (foward) h = lambda * g->nodes[v].dist_to_dest + (1 - lambda) * g->nodes[v].time_to_dest;
                else h = lambda * g->nodes[v].dist_from_source + (1 - lambda) * g->nodes[v].time_from_source;

                double f = new_cost + h;
                if (f > (lambda * dist_best + (1- lambda) * W))
                    continue;

                push(heap, v, f);
            }
        }
        if(u == d || u == s){
            g->nodes[u].active = false;
            continue;
        }

    }

    g->nodes[s].active = true;
    g->nodes[d].active = true;

    if(goal == (uint32_t)-1 && destination != (uint32_t)-1)
        goal = destination;

    AStarResult res;
    uint32_t *path = malloc(sizeof(uint32_t) * n);

    uint32_t length = 0;
    uint32_t v = goal;
    res.dist_cost = 0;
    res.time_cost = 0;
    res.bound_path = 0;

    if(goal != (uint32_t)-1 && cost[goal] != DBL_MAX)
        while (v != start) {
            path[length++] = v;

            if(goal != (uint32_t)-1 && v != start && v != goal)
                if(g->nodes[v].time_from_source + g->nodes[v].time_to_dest > res.bound_path)
                    res.bound_path = g->nodes[v].time_from_source + g->nodes[v].time_to_dest;

            uint32_t u = parent[v];
            if (u == UINT32_MAX) break;

            double add_path_dist = DBL_MAX;
            double add_path_time = DBL_MAX;

            for (uint32_t i = 0; i < g->nodes[u].out_count; i++) {
                Edge e = g->nodes[u].out_edges[i];
                if (e.to == v) {
                    if(add_path_dist > e.distance)
                        add_path_dist = e.distance;
                    if(add_path_time > e.time)
                        add_path_time = e.time;
                }
            }
            res.dist_cost += add_path_dist;
            res.time_cost += add_path_time;

            v = u;
        }

    if(goal != (uint32_t)-1 && cost[goal] == DBL_MAX){
        res.dist_cost = DBL_MAX;
        res.time_cost = DBL_MAX;
    }
    res.path = path;
    res.path_len = length;

    free(cost);
    free(heap->data);
    free(heap);
    free(parent);

    return res;
}
