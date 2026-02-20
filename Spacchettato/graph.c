#include "graph.h"

#include <stdlib.h>
#include <math.h>

#define Pi 0.000000017453292519943295 // Pi / 180 / 1,000,000
#define COEF 61161607.40544 // = Earth radius M * 0.96 * 10

//NOTE: The input coordinates are in nanodegrees (to transform them in degrees they must be multiplied by 1,000,000)
// Length data are given in decimeters (0.1m)

double haversine(int32_t lon1, int32_t lat1, int32_t lon2, int32_t lat2){
    double Delta_lat_squared = pow(fabs(lat1 - lat2) * Pi, 2);
    double cos_mean_lat = cos((lat1 + lat2) * Pi / 2);
    double x = cos_mean_lat * fabs(lon1 - lon2) * Pi;
    return floor(COEF * sqrt(Delta_lat_squared + pow(x, 2)));
}

Graph *createGraph(uint32_t n_nodes, uint32_t n_edges) {
    Graph *g = malloc(sizeof(Graph));
    g->n_nodes = n_nodes;
    g->n_edges = n_edges;
    g->nodes = calloc(n_nodes, sizeof(Node));
    return g;
}

void addEdge(Graph *g, uint32_t from, uint32_t to, double distance, double time) {
    Node *n = &g->nodes[from];
    if (n->out_count == n->out_cap) {
        n->out_cap = n->out_cap ? n->out_cap * 2 : 4;
        n->out_edges = realloc(n->out_edges, n->out_cap * sizeof(Edge));
    }
    n->out_edges[n->out_count++] = (Edge){to, distance, time};

    Node *m = &g->nodes[to];
    if (m->in_count == m->in_cap) {
        m->in_cap = m->in_cap ? m->in_cap * 2 : 4;
        m->in_edges = realloc(m->in_edges, m->in_cap * sizeof(InEdge));
    }
    m->in_edges[m->in_count++] = (InEdge){from, distance, time};
}

void freeGraph(Graph *g) {
    if (!g) return;
    if (g->nodes) {
        for (uint32_t i = 0; i < g->n_nodes; i++) {
            free(g->nodes[i].out_edges);
            free(g->nodes[i].in_edges);
        }
        free(g->nodes);
    }
    free(g);
}
