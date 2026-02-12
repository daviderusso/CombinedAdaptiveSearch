//
// Created by Alessio on 13/04/2025.
//

#ifndef BINARYCOMBINEDHEURISTIC_DIJKSTRA_H
#define BINARYCOMBINEDHEURISTIC_DIJKSTRA_H

#endif //BINARYCOMBINEDHEURISTIC_DIJKSTRA_H

#include "graph.h"
#include "spherical_heuristic.h"
#include "binaryheap.h"

/*
    The dijkstra function (implemented using a binary heap) takes as input:
      - graph: the graph on which to run Dijkstra's algorithm
      - src: source node
      - dest: destination node; if set to -1, the entire tree is computed (i.e., the algorithm runs to completion)
              if dest != -1, the algorithm stops as soon as the destination node is extracted from the queue.
      - dist: array (of size graph->numNodes) in which the minimum distances from src will be stored.
      - pred: array (also of size graph->numNodes) in which the predecessor for each node is stored.
      - useSunCost: if non-zero, the weight used for shortest paths will be arc->cost1; otherwise, arc->cost2.
*/
void dijkstra(Graph *graph, int src, int dest, double *dist, int *pred, int useSunCost) {
    int n = graph->numNodes;
    MinHeap *heap = createMinHeap(n);

    // Initialization of predecessors and distances
    for (int v = 0; v < n; v++) {
        dist[v] = FLT_MAX;
        pred[v] = -1;
        // If the Euclidean distance exceeds the threshold (K*SP_Length), then it should not be inserted into the heap
        // Add a counter for the inserted nodes to adjust the heap size accordingly
        heap->array[v].node = v;
        heap->array[v].dist = FLT_MAX;
        heap->pos[v] = v;
    }

    // Set the distance of the source to zero.
    dist[src] = 0.0f;
    decreaseKey(heap, src, 0.0f);
    heap->size = n;

    while (!isEmpty(heap)) {
        int u = extractMin(heap);

        // If we extracted the destination node, the procedure can end.
        if (dest != -1 && u == dest)
            break;

        Arc *arc = graph->nodes[u].head;
        while (arc != NULL) {
            int v = arc->dest;

            // If useSunCost is enabled, use arc->cost1; otherwise, use arc->cost2.
            double weight;
            if (useSunCost != 0) {
                weight = arc->cost1;
            } else {
                weight = arc->cost2;
            }

            if (isInMinHeap(heap, v) && dist[u] != FLT_MAX && dist[u] + weight < dist[v]) {
                dist[v] = dist[u] + weight;
                pred[v] = u;
                decreaseKey(heap, v, dist[v]);
            }
            arc = arc->next;
        }
    }

    freeMinHeap(heap);
}


/*
    Dijkstra_start function takes as inputs:
      - graph: graph on which to execute Dijkstra
      - src: source node
      - dist: array (of dimension graph->numNodes) in which it will be stored the minimal distances from src.
      - pred: array (of dimension graph->numNodes) in which it will be stored the predecessor of each node.
      - distToDest: array that memorizes a lower bound of the resources needed to reach the destination.
      - W: maximum number of resources available for a feasible path
*/
void dijkstra_start(Graph *graph, int src, double *dist, int *pred, double *distToDest, int W) {
    int n = graph->numNodes;
    double threshold = 0;
    MinHeap *heap = createMinHeap(n);

    // Initialization of predecessors and distances
    for (int v = 0; v < n; v++) {
        // Threshold in case we're exploring the reverse graph
        if (distToDest[src] != 0) {
            threshold = distToDest[v];
        }

        dist[v] = FLT_MAX;
        pred[v] = -1;
        heap->pos[v] = n + 1;

        // Only src and feasible nodes will be inserted in the heap
        if (v == src) {
            dist[v] = 0.0f;
            insertHeap(heap, v, 0.0f);
        } else if (threshold <= W) {
            insertHeap(heap, v, FLT_MAX);
        }
    }

    while (!isEmpty(heap)) {
        // Extract the node whose path length from source is the lowest
        int u = extractMin(heap);

        if (dist[u] != FLT_MAX) {
            Arc *arc = graph->nodes[u].head;
            while (arc != NULL) {
                int v = arc->dest;

                if (isInMinHeap(heap, v)) {
                    // Resources are used as weights
                    double weight = arc->cost2;

                    if (dist[u] + weight < dist[v]) {
                        dist[v] = dist[u] + weight;
                        pred[v] = u;
                        decreaseKey(heap, v, dist[v]);
                    }
                }
                arc = arc->next;
            }
        }
    }

    freeMinHeap(heap);
}

// Reconstructs the path from the source node to the destination node using the predecessor array.
// The path is returned as a dynamic array, and its length is stored in pathLength.
int *reconstructPath(int *pred, int dest, int NumNodes, int *pathLength) {
    int maxNodes = NumNodes;
    int *tempPath = (int *) malloc(maxNodes * sizeof(int));
    int count = 0;
    int current = dest;

    while (current != -1) {
        tempPath[count++] = current;
        current = pred[current];
    }

    // Reverses the path since it was built backwards.
    int *path = (int *) malloc(count * sizeof(int));
    for (int i = 0; i < count; i++) {
        path[i] = tempPath[count - 1 - i];
    }
    free(tempPath);
    *pathLength = count;
    return path;
}

// Computes the sum of costs along a reconstructed path (using the 'cost2' field),
// given the graph, the path, and its length.
double computePathCost2(Graph *graph, int *path, int pathLength) {
    double totalCost = 0.0f;
    // Iterates over consecutive pairs in the path.
    for (int i = 0; i < pathLength - 1; i++) {
        int u = path[i];
        int v = path[i + 1];
        Arc *arc = graph->nodes[u].head;
        double foundCost = -1.0f;
        while (arc != NULL) {
            if (arc->dest == v) {
                foundCost = arc->cost2;
                break;
            }
            arc = arc->next;
        }
        if (foundCost < 0.0f) {
            fprintf(stderr, "Errore: arco da %d a %d non trovato.\n", u, v);
            exit(EXIT_FAILURE);
        }
        totalCost += foundCost;
    }
    return totalCost;
}

// Computes the sum of costs along a reconstructed path (using the 'cost1' field),
// given the graph, the path, and its length.
double computePathCost1(Graph *graph, int *path, int pathLength) {
    double totalCost = 0.0f;
    // Iterates over consecutive pairs in the path.
    for (int i = 0; i < pathLength - 1; i++) {
        int u = path[i];
        int v = path[i + 1];
        Arc *arc = graph->nodes[u].head;
        double foundCost = -1.0f;
        while (arc != NULL) {
            if (arc->dest == v) {
                foundCost = arc->cost1;
                break;
            }
            arc = arc->next;
        }
        if (foundCost < 0.0f) {
            fprintf(stderr, "Errore: arco da %d a %d non trovato.\n", u, v);
            exit(EXIT_FAILURE);
        }
        totalCost += foundCost;
    }
    return totalCost;
}

// Prints the path from the source node to the target node.
void printPath(int source, int target, int *prev) {
    // Computes the length of the path (at most n nodes)
    int pathLength = 0;
    int current = target;

    // First, count how many nodes are in the path
    while (current != -1 && current != source) {
        pathLength++;
        current = prev[current];
    }

    if (current == -1) {
        printf("Nessun cammino da %d a %d.\n", source, target);
        return;
    }

    pathLength++; // include the source node

    // Allocate space for the path
    int *path = (int *) malloc(pathLength * sizeof(int));
    if (!path) {
        perror("Errore allocazione memoria per path");
        exit(EXIT_FAILURE);
    }

    // Reconstruct the path in reverse order
    current = target;
    for (int i = pathLength - 1; i >= 0; i--) {
        path[i] = current;
        current = prev[current];
    }

    // Print in format [node1, node2, ..., noden]
    printf("[");
    for (int i = 0; i < pathLength; i++) {
        printf("%d", path[i]);
        if (i < pathLength - 1)
            printf(", ");
    }
    printf("]\n");

    free(path);
}

void printPathToFile(int source, int target, int *prev, const char *filename) {
    FILE *file = fopen(filename, "w");
    if (!file) {
        perror("Errore apertura file");
        return;
    }

    int pathLength = 0;
    int current = target;

    while (current != -1 && current != source) {
        pathLength++;
        current = prev[current];
    }

    if (current == -1) {
        fprintf(file, "Nessun cammino da %d a %d.\n", source, target);
        fclose(file);
        return;
    }

    pathLength++; // include the source node

    int *path = (int *) malloc(pathLength * sizeof(int));
    if (!path) {
        perror("Errore allocazione memoria per path");
        fclose(file);
        return;
    }

    current = target;
    for (int i = pathLength - 1; i >= 0; i--) {
        path[i] = current;
        current = prev[current];
    }

    // Print to file in the desired format
    fprintf(file, "[");
    for (int i = 0; i < pathLength; i++) {
        fprintf(file, "%d", path[i]);
        if (i < pathLength - 1)
            fprintf(file, ", ");
    }
    fprintf(file, "]\n");

    free(path);
    fclose(file);
}
