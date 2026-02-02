#include <stdio.h>
#include <stdlib.h>
#include "BinaryCombinedHeuristic.h"
#include <string.h>
#include <sys/stat.h>
#include <errno.h>



int main() {
    int TimeLimit = 10;
    const char *filename = "data/W.xy"; // File name of the graph

    clock_t startReadGraph = clock();
    Graph** Grafi = readGraphAndCreateReverse(filename);
    //Graph* graph = readGraphFromFile(filename);
    Graph* graph = Grafi[0];
    Graph* reverseGraph = Grafi[1];
    double TimeReadGraph = (double)(clock() - startReadGraph ) / CLOCKS_PER_SEC;

    // const char *filein = "data/160WInstancesHeur_AP.csv";
    const char *filein = "data/8WInstancesHeur_AP.csv";
    const char *res_dir = "res";
    if (mkdir(res_dir, 0777) != 0) {
        if (errno != EEXIST) {
            perror("Errore creazione cartella res");
            return EXIT_FAILURE;
        }
    }

    char outfile[256];
    char filetime[256];
    snprintf(outfile, sizeof(outfile), "%s/results_heur_W.csv", res_dir);
    snprintf(filetime, sizeof(filetime), "%s/results_times_heur_W.csv", res_dir);


    FILE *fp  = fopen(filein,  "r");
    FILE *fout = fopen(outfile, "w");
    FILE *fouttime = fopen(filetime, "w");
    if (!fp || !fout) {
        perror("Errore apertura file");
        return EXIT_FAILURE;
    }

    char buffer[1024];
    // Skip the header
    if (!fgets(buffer, sizeof(buffer), fp)) {
        fprintf(stderr, "File vuoto o errore in lettura header\n");
        fclose(fp);
        return EXIT_FAILURE;
    }

    // Headers of the CSV
    fprintf(fout, "Source;Target;Budget;runtime_readGraph;runtime_SP_R;runtime_SP_L;Resource_SP_R;Length_SP_R;Length_SP_L;runtime_Red;runtime_BCH;id_time_BCH;Resources_BCH;Length_BCH;BestIter;BestLambda;NumIter\n");
    fprintf(fouttime, "Source;Target;Budget;runtime_readGraph;runtime_SPTreeSource;runtime_SPTreeDest;runtime_Red;runtime_first_lambda;runtime_BCH_total;perc_eliminated_nodes;perc_eliminated_arcs\n");

    char *line = buffer;

    int source, target;
    double W;
    // Read one line at a time
    while (fgets(buffer, sizeof(buffer), fp)) {
        line = buffer;
        // Removes '\n' in queue
        buffer[strcspn(buffer, "\r\n")] = 0;

        // Read data from the line: "%d;%d;%lf"
        if (sscanf(line, "%d;%d;%lf", &source, &target, &W) == 3) {
            printf("\nsource=%d  target=%d  Budget=%lf\n", source, target, W);
        } else {
            fprintf(stderr, "Riga malformata: %s\n", line);
        }

        Results *Risultati = (Results*)malloc(sizeof (Results));

        clock_t startHeur = clock();

        BinaryCombinedHeuristic_AStar(graph, reverseGraph, source, target, W, TimeLimit, Risultati);

        double TimeTotalHeur = (double)(clock() - startHeur ) / CLOCKS_PER_SEC;

        printf("\nIl cammino minimo consuma %lf risorse, e la sua lunghezza è %lf.", Risultati->SP_Length, Risultati->SP_SunLength);
        printf("\nIl cammino Lambda consuma %lf risorse, e la sua lunghezza è %lf. Tempo individuazione: %lf s. Tempo impiegato: %lf s.\n", Risultati->BCH_Length, Risultati->BCH_SunLength, Risultati->idTime_BCH, Risultati->runtime_BCH);

        fprintf(fout, "%d;%d;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%d;%f;%d\n", source,target,W,TimeReadGraph, Risultati->runtime_SP, Risultati->runtime_SPL, Risultati->SP_Length, Risultati->SP_SunLength, Risultati->SPL_Length, Risultati->runtime_red,Risultati->runtime_BCH, Risultati->idTime_BCH, Risultati->BCH_Length, Risultati->BCH_SunLength, Risultati->bestIter, Risultati->bestLambda, Risultati->numIter);
        fflush(fout);

        fprintf(fouttime, "%d;%d;%lf;%lf;%lf;%lf;%lf;%lf;%lf;%f;%f;\n", source,target,W, TimeReadGraph, Risultati->runtime_sptree1, Risultati->runtime_sptree2,Risultati->runtime_red,Risultati->runtime_ALambda, TimeTotalHeur, Risultati->PercRedNodes, Risultati->PercRedArcs);
        fflush(fouttime);

        {
            char detailsFile[256];
            snprintf(detailsFile, sizeof(detailsFile), "%s/%d-%d-%.6f_sol_details.txt", res_dir, source, target, W);
            FILE *fdetails = fopen(detailsFile, "w");
            if (fdetails != NULL) {
                fprintf(fdetails, "index;iter;lambda;timeFound;length;sunLength\n");
                for (int i = 0; i < Risultati->numImprovements; i++) {
                    const Improvement *imp = &Risultati->improvements[i];
                    fprintf(
                        fdetails,
                        "%d;%d;%f;%lf;%lf;%lf\n",
                        i + 1,
                        imp->iter,
                        imp->lambda,
                        imp->timeFound,
                        imp->length,
                        imp->sunLength
                    );
                }
                fclose(fdetails);
            } else {
                perror("Errore apertura file dettagli soluzioni");
            }
        }

        free(Risultati->improvements);
        free(Risultati);


    }

    fclose(fp);
    fclose(fout);


    freeGraph(graph);
    freeGraph(reverseGraph);


    return 0;
}


//VERIFICARE ID - ALESSIO HA AGGIUNTO UNO ALL'id dei nodi. forse si può togliere
