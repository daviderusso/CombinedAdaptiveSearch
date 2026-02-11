//
// Created by Alessio on 18/04/2025.
//

#ifndef BINARYCOMBINEDHEURISTIC_BINARYCOMBINEDHEURISTIC_H
#define BINARYCOMBINEDHEURISTIC_BINARYCOMBINEDHEURISTIC_H

#endif //BINARYCOMBINEDHEURISTIC_BINARYCOMBINEDHEURISTIC_H

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
//#include "Dijkstra.h"
#include "AStar.h"

// Struttura per rappresentare i risultati
typedef struct Improvement {
    double length;
    double sunLength;
    double timeFound;
    float lambda;
    int iter;
} Improvement;

typedef struct Results {
    double runtime_SP;
    double runtime_SPL;
    double SP_sumLength;
    double SP_Length;
    double SPL_Length;
    double runtime_sptree1;
    double runtime_sptree2;
    double runtime_red;
    double runtime_ALambda;
    double runtime_BCH;
    double idTime_BCH;
    double BCH_Length;
    double BCH_SunLength;
    int numIter;
    int bestIter;
    float bestLambda;
    float PercRedNodes;
    float PercRedArcs;
    Improvement *improvements;
    int numImprovements;
    int capImprovements;
} Results;

// HEURISTIC BCH:
// Inputs:
// Graph* graph - Directed graph (in form of adjacency list)
// Graph* reverse_graph - Reverse graph needed to evaluate the shortest path tree from the destination node
// int source, int destination, double W - Instance data: node from which the path starts, node in which it ends and maximum budget for the resources
// Results* Risultati - Structure in which results are stored
// double TimeReadGraph - Time needed to fill the graph and the reverse graph structures.


void BinaryCombinedHeuristic_AStar(Graph *graph, Graph *reverse_graph, int source, int destination, int W,
                                   int TimeLimit, Results *Risultati) {
    // * Preprocessing phase * //

    //Shortest path evaluation (using distances as costs) with A*

    double *distSP = (double *) malloc(graph->numNodes * sizeof(double));
    int *predSP = (int *) malloc(graph->numNodes * sizeof(int));
    clock_t begin_AStarDist = clock();
    astar(graph, source, destination, distSP, predSP, 1); // SP Distances
    clock_t end_AStarDist = clock();
    double time_spent_AstarDist = (double) (end_AStarDist - begin_AStarDist) / CLOCKS_PER_SEC;

    double SP_Cost1_Length = distSP[destination]; //Distance covered by the SP

    Risultati->runtime_SPL = time_spent_AstarDist;

    free(distSP);
    free(predSP);

    //Shortest path evaluation (using resources as costs) with A*

    distSP = (double *) malloc(graph->numNodes * sizeof(double));
    predSP = (int *) malloc(graph->numNodes * sizeof(int));

    clock_t begin_AStarRes = clock();
    astar(graph, source, destination, distSP, predSP, 0); // SP Resources
    clock_t end_AStarRes = clock();
    double time_spent_AstarRes = (double) (end_AStarRes - begin_AStarRes) / CLOCKS_PER_SEC;

    double SP_length = distSP[destination]; //Resources consumed by the SP
    int *pathLength = (int *) malloc(sizeof(int));
    int *SP_path = reconstructPath(predSP, destination, graph->numNodes, pathLength);
    double SP_sunLength = computePathSunCost(graph, SP_path, *pathLength);

    free(SP_path);
    free(pathLength);

    //Graph reduction: for each node v, the minimum resources quantity needed to reach node v from source, and the quantity
    // needed to reach destination from node v are evaluated. If the sum of these two quantities is grater than the budget
    // then node v can be discarded since it cannot be part of a feasible solution (path).
    //These quantities are calculated by solving the shortest path tree problem from the source to each node and then from
    // the destination to each node, using resources as costs.

    clock_t startRed = clock();

    Graph *g_reversed = reverse_graph; //Reverse graph needed for the SPTP from destination to each node

    // Shortest path tree from source to each node
    clock_t startSpTree1 = clock();
    double *distFromSource = (double *) calloc(graph->numNodes, sizeof(double));
    double *distToDest = (double *) calloc(graph->numNodes, sizeof(double));
    int *predFromSource = (int *) malloc(graph->numNodes * sizeof(int));
    dijkstra_start(graph, source, distFromSource, predFromSource, distToDest, W);
    Risultati->runtime_sptree1 = (double) (clock() - startSpTree1) / CLOCKS_PER_SEC;

    // Shortest path tree from destination to each node (on reverse graph)
    clock_t startSpTree2 = clock();
    int *predToDest = (int *) malloc(graph->numNodes * sizeof(int));
    dijkstra_start(g_reversed, destination, distToDest, predToDest, distFromSource, W);
    Risultati->runtime_sptree2 = (double) (clock() - startSpTree2) / CLOCKS_PER_SEC;

    // Free the pointers not needed
    free(distSP);
    free(predSP);
    free(predFromSource);
    free(predToDest);


    // Create reduced graph
    Graph *g_red = reduceGraph(graph, distFromSource, distToDest, W);

    printf("\nGrafo ridotto nodi: %d, archi: %d", g_red->effNodes, g_red->numEdges);
    Risultati->runtime_red = (double) (clock() - startRed) / CLOCKS_PER_SEC;

    // * Start Heuristic * //


    // Initialization Best parameters
    double bestSolLength = SP_length;
    double bestSolSunLength = SP_sunLength;
    float bestLambda = 0;
    double bestTime = 0;
    int num_iter = 0;
    int bestIter = 0;

    Risultati->improvements = NULL;
    Risultati->numImprovements = 0;
    Risultati->capImprovements = 0;


    int check1 = 0;
    bool check = true;

    clock_t startHeur = clock();

    while ((check1 == 0) && ((double) (clock() - startHeur) / CLOCKS_PER_SEC) < TimeLimit) {
        float currLambda = 1;
        num_iter++;

        //Shortest path evaluation (using distances as costs) with A* on the reduced graph

        double *distFromSource_Lambda = (double *) malloc(graph->numNodes * sizeof(double));
        int *predFromSource_Lambda = (int *) malloc(graph->numNodes * sizeof(int));

        clock_t startALambda = clock();

        astar_lambda(g_red, source, destination, distFromSource_Lambda, predFromSource_Lambda, currLambda, distToDest);

        if (check) {
            Risultati->runtime_ALambda = (double) (clock() - startALambda) / CLOCKS_PER_SEC;
            check = false;
        }

        int *lambdaPathLength = (int *) malloc(sizeof(int));
        int *LP_path = reconstructPath(predFromSource_Lambda, destination, g_red->numNodes, lambdaPathLength);
        double LP_length = computePathCost(g_red, LP_path, *lambdaPathLength);
        double LP_sunLength = computePathSunCost(g_red, LP_path, *lambdaPathLength);

        free(distFromSource_Lambda);
        free(predFromSource_Lambda);

        if (LP_length <= W) {
            //Check on feasibility
            if (LP_sunLength < bestSolSunLength || (LP_sunLength == bestSolSunLength && LP_length < bestSolLength)) {
                //Se la soluzione è migliorativa
                //If the path has a lower objective function value than the best solution or it consumes less resources, then
                // it becomes the new incumbent solution
                bestSolSunLength = LP_sunLength;
                bestSolLength = LP_length;
                bestIter = num_iter;
                bestLambda = currLambda;
                bestTime = (double) (clock() - startHeur) / CLOCKS_PER_SEC;
                // printf("\nSoluzione trovata di costo %lf", bestSolLength);
                // printf("\nSoluzione trovata di costo sun %lf", bestSolSunLength);
                if (Risultati->numImprovements == Risultati->capImprovements) {
                    int newCap = (Risultati->capImprovements == 0) ? 4 : (Risultati->capImprovements * 2);
                    Improvement *newArr = (Improvement *) realloc(Risultati->improvements,
                                                                  (size_t) newCap * sizeof(Improvement));
                    if (newArr != NULL) {
                        Risultati->improvements = newArr;
                        Risultati->capImprovements = newCap;
                    }
                }
                if (Risultati->numImprovements < Risultati->capImprovements) {
                    Improvement *slot = &Risultati->improvements[Risultati->numImprovements];
                    slot->length = bestSolLength;
                    slot->sunLength = bestSolSunLength;
                    slot->timeFound = bestTime;
                    slot->lambda = bestLambda;
                    slot->iter = bestIter;
                    Risultati->numImprovements++;
                }
            }
            check1 = 1; //At least a feasible solution was found
        } else {
            //If the path is not feasible, a binary search for the best value of lambda starts

            // In order to further reduce the graph and to obtain a different solution, the farthest node (the node which consumes
            // the most resources in the path) is sought, and the mininum quantity of resources it needs is stored in newW.
            float newW = maxSumDistanceOnPath(LP_path, *lambdaPathLength, distFromSource, distToDest);

            // Lambda parameters initialization
            float lambda_sup = 1;
            currLambda = 0.5f;
            float lambda_inf = 0;
            int maxIter = 10;
            int curr_iter = 0;

            while ((curr_iter < maxIter) && ((double) (clock() - startHeur) / CLOCKS_PER_SEC < TimeLimit)) {
                curr_iter++;
                num_iter++;

                //Shortest path evaluation (using a convex combination of distances and resources as costs) with A*

                distFromSource_Lambda = (double *) malloc(graph->numNodes * sizeof(double));
                predFromSource_Lambda = (int *) malloc(graph->numNodes * sizeof(int));
                astar_lambda(g_red, source, destination, distFromSource_Lambda, predFromSource_Lambda, currLambda,
                             distToDest);

                free(LP_path);
                free(lambdaPathLength);

                lambdaPathLength = (int *) malloc(sizeof(int));
                LP_path = reconstructPath(predFromSource_Lambda, destination, g_red->numNodes, lambdaPathLength);
                LP_length = computePathCost(g_red, LP_path, *lambdaPathLength);
                LP_sunLength = computePathSunCost(g_red, LP_path, *lambdaPathLength);

                free(distFromSource_Lambda);
                free(predFromSource_Lambda);

                if (LP_length <= W) {
                    // Check on feasibility
                    if (LP_sunLength < bestSolSunLength || (
                            LP_sunLength == bestSolSunLength && LP_length < bestSolLength)) {
                        //Se la soluzione è migliorativa
                        //If the path has a lower objective function value than the best solution or it consumes less resources, then
                        // it becomes the new incumbent solution
                        bestSolSunLength = LP_sunLength;
                        bestSolLength = LP_length;
                        bestIter = num_iter;
                        bestLambda = currLambda;
                        bestTime = (double) (clock() - startHeur) / CLOCKS_PER_SEC;
                        // printf("\nSoluzione trovata di costo %lf", bestSolLength);
                        // printf("\nSoluzione trovata di costo sun %lf", bestSolSunLength);
                        if (Risultati->numImprovements == Risultati->capImprovements) {
                            int newCap = (Risultati->capImprovements == 0) ? 4 : (Risultati->capImprovements * 2);
                            Improvement *newArr = (Improvement *) realloc(
                                Risultati->improvements, (size_t) newCap * sizeof(Improvement));
                            if (newArr != NULL) {
                                Risultati->improvements = newArr;
                                Risultati->capImprovements = newCap;
                            }
                        }
                        if (Risultati->numImprovements < Risultati->capImprovements) {
                            Improvement *slot = &Risultati->improvements[Risultati->numImprovements];
                            slot->length = bestSolLength;
                            slot->sunLength = bestSolSunLength;
                            slot->timeFound = bestTime;
                            slot->lambda = bestLambda;
                            slot->iter = bestIter;
                            Risultati->numImprovements++;
                        }
                    }
                    // Increase the value of lambda in order to improve the quality of the solution
                    lambda_inf = currLambda;
                    currLambda = (lambda_sup - lambda_inf) / 2 + lambda_inf;
                } else {
                    //If path is infeasible

                    // Decrease the value of lambda in order to restore feasibility
                    lambda_sup = currLambda;
                    currLambda = (lambda_sup - lambda_inf) / 2 + lambda_inf;
                }

                //printPath(source, destination, predFromSource_Lambda);
                //printPathToFile(source, destination, predFromSource_Lambda, "camminoC.txt");
            }

            // * Graph Reduction*

            free(LP_path);
            free(lambdaPathLength);

            if (((double) (clock() - startHeur) / CLOCKS_PER_SEC) < TimeLimit) {
                // Exclude every node v which has a minimum resource quantity needed greater than
                Graph *reduced = reduceGraph(g_red, distFromSource, distToDest, newW);

                // Free the old graph since it's not needed.
                freeGraph(g_red);

                // Overwrite the original pointer.
                g_red = reduced;
            }

            //printf("\nNumero archi dopo riduzione:%d", g_red->numEdges);
        }
    }

    clock_t endHeur = clock();

    Risultati->PercRedNodes = (float) (graph->effNodes - g_red->effNodes) / (float) (graph->effNodes) * 100;
    Risultati->PercRedArcs = (float) (graph->numEdges - g_red->numEdges) / (float) (graph->numEdges) * 100;
    //printf("\nNumero archi a fine euristica:%d", g_red->numEdges);

    freeGraph(g_red);
    free(distFromSource);
    free(distToDest);

    // Save results in the structure
    Risultati->runtime_SP = time_spent_AstarRes;
    Risultati->SP_sumLength = SP_length; //calcolo sp1 prendo distanza (min distanza)
    Risultati->SP_Length = SP_sunLength; //calcolo sp1 prendo risorsa (min distanza)
    Risultati->SPL_Length = SP_Cost1_Length; //calcolo sp2 prendo distanza (min risorsa)
    Risultati->runtime_BCH = (double) (endHeur - startHeur) / CLOCKS_PER_SEC;
    Risultati->idTime_BCH = bestTime;
    Risultati->BCH_Length = bestSolLength;
    Risultati->BCH_SunLength = bestSolSunLength;
    Risultati->bestLambda = bestLambda;
    Risultati->bestIter = bestIter;
    Risultati->numIter = num_iter;
}
