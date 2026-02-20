#ifndef GRAPH_H
#define GRAPH_H

#include <stdint.h>
#include <stdbool.h>

typedef struct {
    uint32_t to;
    double distance;
    double time;
} Edge;

typedef struct {
    uint32_t from;
    double distance;
    double time;
} InEdge;

typedef struct {
    Edge *out_edges;
    InEdge *in_edges;
    uint32_t out_count, in_count;
    uint32_t out_cap, in_cap;
    bool active;
    double x;
    double y;
    double dist_from_source;
    double dist_to_dest;
    double time_from_source;
    double time_to_dest;
} Node;

typedef struct {
    Node *nodes;
    uint32_t n_nodes;
    uint32_t n_edges;
} Graph;

double haversine(int32_t lon1, int32_t lat1, int32_t lon2, int32_t lat2);
Graph *createGraph(uint32_t n_nodes, uint32_t n_edges);
void addEdge(Graph *g, uint32_t from, uint32_t to, double distance, double time);
void freeGraph(Graph *g);

#endif
