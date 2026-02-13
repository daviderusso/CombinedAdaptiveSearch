//
// Created by Alessio on 12/04/2025.
//

#ifndef BINARYCOMBINEDHEURISTIC_GRAPH_H
#define BINARYCOMBINEDHEURISTIC_GRAPH_H

#endif //BINARYCOMBINEDHEURISTIC_GRAPH_H

#include <stdlib.h>
#include <float.h>

// I represent the graph as an adjacency list: for each node, I have a list of all nodes adjacent to it.
// Each element of this list represents an edge with its properties.

// NOTE: Node indexing for our heuristic starts from 0, while for the exact algorithm it starts from 1!

// Represents an edge in the adjacency list.
typedef struct Arc {
    int dest; // Destination node.
    double cost1; // First cost (cost to minimize).
    double cost2; // Second cost (cost constrained by a bound).
    struct Arc *next; // Pointer to the next edge.
} Arc;

// Graph node: contains the head of the list of outgoing edges.
typedef struct {
    Arc *head;
    int x, y; // Node coordinates.
} Node;

// Graph structure: number of nodes, number of edges, number of usable nodes, and dynamic array of nodes.
typedef struct {
    int numNodes;
    int numEdges;
    int effNodes;
    Node *nodes;
} Graph;

// Creates a graph with numNodes nodes.
Graph *createGraph(int numNodes) {
    Graph *g = (Graph *) malloc(sizeof(Graph));
    if (g == NULL) {
        fprintf(stderr, "Errore nell'allocazione del grafo\n");
        exit(EXIT_FAILURE);
    }
    g->numNodes = numNodes;
    g->effNodes = numNodes;
    g->nodes = (Node *) calloc(numNodes, sizeof(Node));
    if (g->nodes == NULL) {
        fprintf(stderr, "Errore nell'allocazione dei nodi\n");
        exit(EXIT_FAILURE);
    }
    return g;
}

// Adds an edge to the graph (inserting it at the head of the adjacency list of node u).
void addArc(Graph *g, int u, int v, double cost1, double cost2) {
    Arc *newArc = (Arc *) malloc(sizeof(Arc));
    if (newArc == NULL) {
        fprintf(stderr, "Errore nell'allocazione di un arco\n");
        exit(EXIT_FAILURE);
    }
    newArc->dest = v;
    newArc->cost1 = cost1;
    newArc->cost2 = cost2;
    newArc->next = g->nodes[u].head;
    g->nodes[u].head = newArc;
}


/*
    The function maxSumDistanceOnPath computes, for each node in the path,
    the sum: d(source, v) + d(v, destination)
    and returns the maximum value obtained.

    Parameters:
    - path: an array of integers containing the node indices along the path.
    - pathLength: length of the path array (number of nodes in the path).
    - distSource: array of doubles containing the shortest distances from the source to each node.
    - distDest: array of doubles containing the shortest distances from the destination to each node.

    The function assumes the distances were computed using Dijkstra or A* and are ≥ 0.
*/
double maxSumDistanceOnPath(int *path, int pathLength, double *distSource, double *distDest) {
    double maxSum = 0.0f; // Dato che le distanze sono ≥ 0, inizializziamo a zero.

    for (int i = 0; i < pathLength; i++) {
        int node = path[i];
        double currentSum = distSource[node] + distDest[node];
        if (currentSum > maxSum) {
            maxSum = currentSum;
        }
    }
    return maxSum;
}

/*
    The function reduceGraph creates a new graph (with the same node numbering) in which:
      - Only "useful" nodes are kept, i.e., those for which:
            d(source,v) + d(v,destination) <= K.
      - For each edge u -> v, the edge is inserted only if both nodes satisfy the condition.
    The d_source and d_dest arrays must have been previously computed (with Dijkstra or A* on the original and reverse graphs).
*/
Graph *reduceGraph(Graph *original, double *d_source, double *d_dest, double W) {
    int n = original->numNodes;
    Graph *reduced = createGraph(n);
    int edgeCount = 0;
    reduced->effNodes = 0; // Reset the effective number of nodes used by the reduced graph.

    for (int u = 0; u < n; u++) {
        // Check if node u is feasible.
        if (d_source[u] > W || d_dest[u] > W || d_source[u] + d_dest[u] > W)
            continue; // Skip adding the node if it doesn't meet the criteria.

        reduced->nodes[u].x = original->nodes[u].x;
        reduced->nodes[u].y = original->nodes[u].y;

        Arc *arc = original->nodes[u].head;
        reduced->effNodes = reduced->effNodes + 1;
        while (arc) {
            int v = arc->dest;
            // Insert the edge only if the adjacent node v also satisfies the condition.
            if (d_source[v] <= W && d_dest[v] <= W && d_source[v] + d_dest[v] <= W && d_source[u] + arc->cost2 + d_dest[
                    v] <= W) {
                addArc(reduced, u, v, arc->cost1, arc->cost2);
                edgeCount++;
            }
            arc = arc->next;
        }
    }

    reduced->numEdges = edgeCount;
    return reduced;
}

void freeGraph(Graph *g) {
    if (g == NULL)
        return;

    // Frees every adjacency list.
    for (int i = 0; i < g->numNodes; i++) {
        Arc *current = g->nodes[i].head;
        while (current != NULL) {
            Arc *next = current->next;
            free(current);
            current = next;
        }
    }
    // Frees the node array and the graph structure.
    free(g->nodes);
    free(g);
}

// Reads a graph and simultaneously creates its reverse graph.
Graph **readGraphAndCreateReverse(const char *filename) {
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) {
        perror("Errore nell'apertura del file");
        exit(EXIT_FAILURE);
    }

    int numNodes, numEdges;

    // Reads the first line which contains "nodes" and "edges"
    if (fscanf(fp, "nodes %d edges %d\n", &numNodes, &numEdges) != 2) {
        fprintf(stderr, "Formato file errato nella prima linea\n");
        exit(EXIT_FAILURE);
    }
    Graph **grafi = malloc(2 * sizeof(Graph *));

    Graph *g = createGraph(numNodes);
    Graph *g_rev = createGraph(numNodes);

    grafi[0] = g;
    grafi[1] = g_rev;

    g->numEdges = numEdges;
    g_rev->numEdges = numEdges;

    // Reads the node lines
    // Each line in the file has the format:
    //   v <index> <x_coordinate> <y_coordinate>
    char type;
    int index;
    int x, y;
    for (int i = 0; i < numNodes; i++) {
        if (fscanf(fp, " %c", &type) != 1) {
            fprintf(stderr, "Errore nella lettura del tipo di riga\n");
            exit(EXIT_FAILURE);
        }
        if (type != 'v') {
            fprintf(stderr, "Errore: atteso 'v', ottenuto '%c'\n", type);
            exit(EXIT_FAILURE);
        }
        if (fscanf(fp, " %d %d %d\n", &index, &x, &y) != 3) {
            fprintf(stderr, "Formato file errato nella riga dei nodi\n");
            exit(EXIT_FAILURE);
        }
        g->nodes[index].x = x;
        g->nodes[index].y = y;
        g_rev->nodes[index].x = x;
        g_rev->nodes[index].y = y;
    }

    // Reads the edge lines
    // Each edge line has the format:
    //   e <source_node> <destination_node> <cost1> <cost2>
    int u, v, cost1, cost2;
    for (int i = 0; i < numEdges; i++) {
        if (fscanf(fp, " %c", &type) != 1) {
            fprintf(stderr, "Errore nella lettura del tipo di riga per un arco\n");
            exit(EXIT_FAILURE);
        }
        if (type != 'e') {
            fprintf(stderr, "Errore: atteso 'e', ottenuto '%c'\n", type);
            exit(EXIT_FAILURE);
        }
        if (fscanf(fp, " %d %d %d %d\n", &u, &v, &cost1, &cost2) != 4) {
            fprintf(stderr, "Formato file errato nella riga degli archi\n");
            exit(EXIT_FAILURE);
        }

        if (cost1 < 0 || cost2 < 0) {
            printf("ERRORE - costi archi con valori negativi - C1: %lf - C2: %lf\n", cost1, cost2);
        }

        addArc(g, u, v, cost1, cost2);
        addArc(g_rev, v, u, cost1, cost2);
    }

    fclose(fp);
    return grafi;
}
