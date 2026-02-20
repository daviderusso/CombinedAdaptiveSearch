//
// Created by Alessio on 24/04/2025.
//

#ifndef BINARYCOMBINEDHEURISTIC_ASTAR_H
#define BINARYCOMBINEDHEURISTIC_ASTAR_H


//#include "graph.h"
//#include "spherical_heuristic.h"
//#include "binaryheap.h"  //Tutti inclusi in Dijkstra
#include "Dijkstra.h"
#include <stdbool.h>

//
// A* functions
//

// This function takes as input: the graph g, the source src, the destination dest, and two arrays of length n (number of nodes)
// which will respectively contain the distance from the source to node i (dist[i]) and the predecessor of i (pred[i]).

void astar(Graph *g, int src, int dest, double *dist, int *pred, int useCost1) {
    int n = g->numNodes;
    MinHeap *open = createMinHeap(n);
    bool *inOpen = calloc(n, sizeof(bool));

    // Initialize dist and pred
    for (int v = 0; v < n; v++) {
        dist[v] = DBL_MAX;
        pred[v] = -1;
    }

    // Insert only src into the open-set
    dist[src] = 0.0;
    double h0 = 0;
    insertHeap(open, src, h0);
    inOpen[src] = true;

    while (!isEmpty(open)) {
        int u = extractMin(open);
        inOpen[u] = false;
        if (u == dest) break;

        for (Arc *a = g->nodes[u].head; a; a = a->next) {
            int v = a->dest;
            double w;
            if (useCost1 != 0) {
                w = a->cost1;
            } else {
                w = a->cost2;
            }
            double tentative_g = dist[u] + w;

            if (tentative_g < dist[v]) {
                dist[v] = tentative_g;
                pred[v] = u;
                double f = tentative_g + 0;

                if (inOpen[v]) {
                    decreaseKey(open, v, f);
                } else {
                    insertHeap(open, v, f);
                    inOpen[v] = true;
                }
            }
        }
    }

    free(inOpen);
    freeMinHeap(open);
}

void astar_lambda(Graph *g, int src, int dest, double *dist, int *pred, double lambda, double *distToDest) {
    int n = g->numNodes;
    MinHeap *open = createMinHeap(n);
    bool *inOpen = calloc(n, sizeof(bool));

    // Initialize dist and pred
    for (int v = 0; v < n; v++) {
        dist[v] = DBL_MAX;
        pred[v] = -1;
    }

    // Insert only src into the open-set
    dist[src] = 0.0;
    // double h0 = (lambda)*spherical_heuristic_get_explicit(g, src, dest);
    double h0 = 0;
    insertHeap(open, src, h0);
    inOpen[src] = true;

    while (!isEmpty(open)) {
        int u = extractMin(open);
        inOpen[u] = false;
        if (u == dest) break;

        for (Arc *a = g->nodes[u].head; a; a = a->next) {
            int v = a->dest;
            double w;
            // Cost selection based on the lambda value
            w = (lambda) * a->cost1 + (1 - lambda) * a->cost2;

            double tentative_g = dist[u] + w;

            if (tentative_g < dist[v]) {
                dist[v] = tentative_g;
                pred[v] = u;
                double f = tentative_g + (lambda) * spherical_heuristic_get_explicit(g, v, dest) + (1 - lambda) *
                           distToDest[v];
                // printf("\nDistance = %lf", spherical_heuristic_get_explicit(g, v, dest));
                if (inOpen[v]) {
                    decreaseKey(open, v, f);
                } else {
                    insertHeap(open, v, f);
                    inOpen[v] = true;
                }
            }
        }
    }

    free(inOpen);
    freeMinHeap(open);
}

#endif //BINARYCOMBINEDHEURISTIC_ASTAR_H
