#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
#include <stdint.h>
#include <float.h>
#include <inttypes.h>
#include <sys/stat.h>
#include <errno.h>

#include "astar.h"
#include "graph.h"
#include "improvements.h"

#define INF INT_MAX

#define BUF_SIZE (1 << 20)  // 1 MB

// ---------- Parser ----------
void parse_args(int argc, char *argv[],
                char **input_path, uint32_t *s, uint32_t *d,
                uint64_t *W, double *time_limit, int *run_reduction_heuristic, int *max_iterations, int *runall,
                char **inputrunall, char **outdir)
{
    *input_path = NULL;
    *s = *d = 0;
    *W = 0;
    *time_limit = 0.0;
    *run_reduction_heuristic = 1;
    *max_iterations = 0;   
    *runall = 0;
    *inputrunall = NULL;
    *outdir = NULL;
    for (int i = 1; i < argc; i++) {
        if (!strcmp(argv[i], "--input")) {
            *input_path = argv[++i];
        }
        else if (!strcmp(argv[i], "--s")) {
            *s = (uint32_t) atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "--d")) {
            *d = (uint32_t) atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "--W")) {
            *W = strtoull(argv[++i], NULL, 10);
        }
        else if (!strcmp(argv[i], "--tl")) {
            *time_limit = atof(argv[++i]);
        }
        else if (!strcmp(argv[i], "--redh")) {
            *run_reduction_heuristic = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "--nit")) {
            *max_iterations = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "--runall")) {
            *runall = atoi(argv[++i]);
        }
        else if (!strcmp(argv[i], "--inputrunall")) {
            *inputrunall = argv[++i];
        }
        else if (!strcmp(argv[i], "--outdir")) {
            *outdir = argv[++i];
        }
        else {
            fprintf(stderr, "Unknown argument: %s\n", argv[i]);
            exit(EXIT_FAILURE);
        }
    }

    if (!*input_path) {
        fprintf(stderr,
            "Usage: ./compute_paths "
            "--input file --s src --d dest --W bound "
            "[--tl time_limit_seconds] [--nit max_iterations] [--runall run all instances] "
            "[--inputrunall file run all instances] [--outdir output_directory]\n"
        );
        exit(EXIT_FAILURE);
    }
}

// ---------- Main ----------
int main(int argc, char *argv[]) {
    clock_t start_read = clock();
    clock_t end_read = clock();
    char *filename;
    char *inputrunall;
    char *outdir;
    uint32_t s, d;
    double lambda = 0.5;
    uint64_t W_read;
    double time_limit;
    int run_reduction_heuristic = 1;
    int max_iterations;
    int runall;
    double dist_spt, time_spt, dist_spd, time_spd;
    
    parse_args(argc, argv, &filename, &s, &d, &W_read, &time_limit, &run_reduction_heuristic,
               &max_iterations, &runall, &inputrunall, &outdir);

    if (!outdir) {
        fprintf(stderr, "Error: --outdir is required\n");
        return 1;
    }
    if (mkdir(outdir, 0755) != 0 && errno != EEXIST) {
        perror("Error creating output directory");
        return 1;
    }

    start_read = clock();

    FILE *f = fopen(filename, "r");
    if (!f) { perror("Error opening file"); return 1; }

    uint32_t n_nodes, n_edges;
    fscanf(f, "nodes %u edges %u\n", &n_nodes, &n_edges);
    fflush(f);
    fclose(f);

    Graph *g = createGraph(n_nodes, n_edges);

    FILE *fp = fopen(filename, "rb");
    if (!fp) {
        perror("fopen");
        return 1;
    }

    static unsigned char buf[BUF_SIZE];
    size_t n;



    int idx = 0;
    uint32_t id = 0;
    uint32_t node_i = 0;
    uint32_t edge_i = 0;
    uint32_t from;
    uint32_t to;
    uint32_t dist;
    uint32_t time;




    int32_t  num = 0;
    int32_t  vals[4];
    char header[256];

    int      sign = 1;
    int      reading = 0;
    int      is_node = 1;

    fgets(header, sizeof(header), fp);
    sscanf(header, "nodes %u edges %u", &n_nodes, &n_edges);

    while ((n = fread(buf, 1, BUF_SIZE, fp)) > 0) {
        for (size_t i = 0; i < n; i++) {
            unsigned char c = buf[i];

            if (c == 'v' || c == 'e') {
                idx = 0;
                continue;
            }

            if (c == '-') {
                sign = -1;
                reading = 1;
                continue;
            }

            if (c >= '0' && c <= '9') {
                num = num * 10 + (c - '0');
                reading = 1;
                continue;
            }

            if (reading) {
                vals[idx++] = sign * num;
                num = 0;
                sign = 1;
                reading = 0;

                if (is_node && idx == 3) {
                    id = (uint32_t)vals[0];
                    g->nodes[id].x  = vals[1];
                    g->nodes[id].y  = vals[2];
		    g->nodes[id].active = true;
		    g->nodes[id].time_from_source = 0;
		    g->nodes[id].time_to_dest = 0;
		    g->nodes[id].dist_from_source = 0; //haversine(x[s], y[s], x[id], y[id]);
		    g->nodes[id].dist_to_dest = 0; //haversine(x[d], y[d], x[id], y[id]);
		    //printf("x = %" PRId32 " y = %" PRId32 "\n", x[id], y[id]);
		    //getchar();
                    node_i++;
                    idx = 0;
                    if (node_i == n_nodes) is_node = 0;
                }
                else if (!is_node && idx == 4) {
		  from = (uint32_t)vals[0];
		  to = (uint32_t)vals[1];
		  dist = (uint32_t)vals[2];
		  time = (uint32_t)vals[3];
		  //printf("from = %" PRIu32 " to = %" PRIu32 " dist = %" PRIu32 " time = %" PRIu32 "\n", from, to, dist, time);
		  //getchar();

		  addEdge(g, from, to, (double) dist, (double) time);
		  if(edge_i < n_nodes){
		    g->nodes[edge_i].dist_from_source = haversine(g->nodes[s].x, g->nodes[s].y, g->nodes[edge_i].x, g->nodes[edge_i].y);
		    g->nodes[edge_i].dist_to_dest = haversine(g->nodes[d].x, g->nodes[d].y, g->nodes[edge_i].x, g->nodes[edge_i].y);
		  }
		  edge_i++;
		  idx = 0;
                }
            }
        }
    }


    if (reading) {
        vals[idx++] = sign * num;

        if (is_node && idx == 3) {
	  id = vals[0];
	  g->nodes[id].x  = (int32_t)vals[1];
	  g->nodes[id].y  = (int32_t)vals[2];
	  g->nodes[id].active = true;
	  g->nodes[id].time_from_source = 0;
	  g->nodes[id].time_to_dest = 0;
	  g->nodes[id].dist_from_source = 0; //haversine(x[s], y[s], x[id], y[id]);
	  g->nodes[id].dist_to_dest = 0; //haversine(x[d], y[d], x[id], y[id]);
	  node_i++;
        } else if (!is_node && idx == 4) {
	  from = (uint32_t)vals[0];
	  to   = (uint32_t)vals[1];
	  dist = (uint32_t)vals[2];
	  time = (uint32_t)vals[3];
	  addEdge(g, from, to, (double) dist, (double) time);
	  if(edge_i < n_nodes){
	    g->nodes[edge_i].dist_from_source = haversine(g->nodes[s].x, g->nodes[s].y, g->nodes[edge_i].x, g->nodes[edge_i].y);
	    g->nodes[edge_i].dist_to_dest = haversine(g->nodes[d].x, g->nodes[d].y, g->nodes[edge_i].x, g->nodes[edge_i].y);
	  }
	  edge_i++;
        }

    }

    fclose(fp);

    for (uint32_t i = n_edges; i < n_nodes; i++) {
      g->nodes[i].dist_from_source = haversine(g->nodes[s].x, g->nodes[s].y, g->nodes[i].x, g->nodes[i].y);
      g->nodes[i].dist_to_dest = haversine(g->nodes[d].x, g->nodes[d].y, g->nodes[i].x, g->nodes[i].y);
    }

    //printf("Graph loaded: %u nodes.\n", n_nodes);

    end_read = clock();
    //printf("Time read graph: %lf\n", ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)));


    FILE *input_runall;
    int num_ist = 1;
    if(runall == 1){
      input_runall = fopen(inputrunall,"r");
      if (!input_runall) { perror("Error opening file"); return 1; }
      fscanf(input_runall, "%d\n", &num_ist);
    }

    char results_path[512];
    if (snprintf(results_path, sizeof(results_path), "%s/%s", outdir, "results.txt") >= (int)sizeof(results_path)) {
        fprintf(stderr, "Output path too long for results.txt\n");
        return 1;
    }
    FILE *fout = fopen(results_path,"w");
    fprintf(fout, "Source Target Budget runtime_readGraph runtime_SP Length_SP SunLength_SP Length_SPL runtime_Red runtime_BCH id_time_BCH Length_BCH SunLength_BCH BestIter BestLambda NumIter total_runtime total_runtime_best_solution\n");


    for( int ist = 0; ist < num_ist; ist++){

      printf("Istance %d source %d target %d budget %lu\n\n", ist+1, s, d, W_read);

      if(runall == 1)
	fscanf(input_runall, "%u %u %lu\n", &s, &d, &W_read);

      char improvements_filename[512];
      if (snprintf(improvements_filename, sizeof(improvements_filename),
                   "%s/%" PRIu32 "-%" PRIu32 "-%" PRIu64 "_sol_details.txt",
                   outdir, s, d, (uint64_t) W_read) >= (int)sizeof(improvements_filename)) {
          fprintf(stderr, "Output path too long for improvements file\n");
          return 1;
      }
      ImprovementList improvements;
      improvements_init(&improvements);

      double W = (double) W_read;

	for (uint32_t i = 0; i < n_nodes; i++) {
	  g->nodes[i].active = true;
	  g->nodes[i].time_from_source = 0;
	  g->nodes[i].time_to_dest = 0;
	  g->nodes[i].dist_from_source = haversine(g->nodes[s].x, g->nodes[s].y, g->nodes[i].x, g->nodes[i].y);
	  g->nodes[i].dist_to_dest = haversine(g->nodes[d].x, g->nodes[d].y, g->nodes[i].x, g->nodes[i].y);
	}

      clock_t start_exec = clock();
      clock_t end_exec = clock();

      AStarResult path_inf = {0};
      AStarResult path_best = {0};
      uint32_t remaining = 0;

      clock_t start_time_red = clock();
      clock_t end_time_red = clock();


      //printf("Compute shortest-path distance...\n");
      clock_t start_time_spd = clock();
      clock_t end_time_spd = clock();


      AStarResult astar_spd = a_star_with_bound(g, s, d, DBL_MAX, DBL_MAX, true, s, d, W, DBL_MAX, 1, d); // forward, time


      end_time_spd = clock();
      dist_spd = astar_spd.dist_cost;
      time_spd = astar_spd.time_cost;
      if(astar_spd.time_cost > W){
	//printf("Shortest-path distance is infeasible.\n");
	//printf("Cost (distance): %lf, Cost (time): %lf\n", astar_spd.dist_cost, astar_spd.time_cost);

	free_astar_result(&path_inf);
	path_inf = copy_res(&astar_spd);

	if(path_inf.time_cost == DBL_MAX){
	  //printf("Exit: Graph is no more connected.\n");
	  fprintf(fout, "%u %u %.0lf %.3f 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0 0 0 1 0 %.3f %.3f\n", s, d, (double)W_read, ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)));
	  improvements_write_file(improvements_filename, &improvements);
	  improvements_free(&improvements);
	  continue;
	}

      }else{
	//printf("Shortest-path distance is feasible.\n");
	free_astar_result(&path_best);
	path_best = copy_res(&astar_spd);
	improvements_push(&improvements, path_best.dist_cost, path_best.time_cost,
			  elapsed_seconds(start_exec), 1.0, 0);

	end_exec = clock();
	//printf("Best solution found.\n");
	//printf("Cost (distance): %lf, Cost (time): %lf best_iter_red: 0 best_iter_lambda 0 time read graph: %.3f time-preprocessing: %.3f time heuristic 0.0\n", path_best.dist_cost, path_best.time_cost, ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), ((double)((double)(end_exec- start_exec) / CLOCKS_PER_SEC)));
	fprintf(fout, "%u %u %.0lf %.3f 0.0 0.0 0.0 %.0f 0.0 0.0 0.0 0 0 0 1 0 %.3f %.3f\n", s, d, (double)W_read, ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), dist_spd, ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)));
	improvements_write_file(improvements_filename, &improvements);
	improvements_free(&improvements);
	continue;

      }
      free_astar_result(&astar_spd);

      //printf("\nCompute shortest-path resource...\n");

      clock_t start_time_spt = clock();
      clock_t end_time_spt = clock();

      //AStarResult astar_spr = a_star_with_bound(g, s, -1, DBL_MAX, DBL_MAX, true, s, d, W, DBL_MAX, 0, d); // forward, time
      AStarResult astar_spr = a_star_with_bound(g, s, -1, W, DBL_MAX, true, s, d, W, DBL_MAX, 0, d); // forward, time


      end_time_spt = clock();
      dist_spt = astar_spr.dist_cost;
      time_spt = astar_spr.time_cost;

      if(astar_spr.time_cost <= W){

	//printf("Shortest-path time is feasible.\n");
	//printf("Cost (distance): %lf, Cost (time): %lf\n", astar_spr.dist_cost, astar_spr.time_cost);
	free_astar_result(&path_best);
	path_best = copy_res(&astar_spr);
	improvements_push(&improvements, path_best.dist_cost, path_best.time_cost,
			  elapsed_seconds(start_exec), 0.0, 0);

      }else if(astar_spr.time_cost > W){
	//printf("Exit: Shortest-path time is infeasible.\n");
	fprintf(fout, "%u %u %.0lf %.3f 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0 0 0 1 0 %.3f %.3f\n", s, d, (double)W_read, ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)));
	improvements_write_file(improvements_filename, &improvements);
	improvements_free(&improvements);
	continue;
      }else{
	//printf("Exit: Graph is no more connected.\n");
	fprintf(fout, "%u %u %.0lf %.3f 0.0 0.0 0.0 0.0 0.0 0.0 0.0 0 0 0 1 0 %.3f %.3f\n", s, d, (double)W_read, ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)));
	improvements_write_file(improvements_filename, &improvements);
	improvements_free(&improvements);
	continue;

      }
      free_astar_result(&astar_spr);

      // Conta i nodi rimasti

      remaining = 0;
      for (uint32_t i = 0; i < g->n_nodes; i++)
	if(g->nodes[i].active)
	if(i != s && i != d){
	  if (g->nodes[i].time_to_dest + g->nodes[i].time_from_source > 0 && g->nodes[i].time_to_dest + g->nodes[i].time_from_source <=  W && g->nodes[i].dist_to_dest + g->nodes[i].dist_from_source < path_best.dist_cost) remaining++;
	  else g->nodes[i].active = false;
	  if (g->nodes[i].time_from_source == 0) g->nodes[i].active = false;

	}
      //printf("Graph reduced: %u nodes remain active.\n", remaining);

      //printf("\nReducing graph backward using A* (time) with budget W = %lf...\n", W);

      //AStarResult astar_backward_reduction = a_star_with_bound(g, d, -1, DBL_MAX, DBL_MAX, false, s, d, W, path_best.dist_cost, 0, d); // forward, time
      AStarResult astar_backward_reduction = a_star_with_bound(g, d, -1, W, path_best.dist_cost, false, s, d, W, path_best.dist_cost, 0, d); // forward, time

      free_astar_result(&astar_backward_reduction);

      // Conta i nodi rimasti

      remaining = 0;
      for (uint32_t i = 0; i < g->n_nodes; i++)
	if(g->nodes[i].active)
	if(i != s && i != d){
	  if (g->nodes[i].time_to_dest + g->nodes[i].time_from_source > 0 && g->nodes[i].time_to_dest + g->nodes[i].time_from_source <=  W && g->nodes[i].dist_to_dest + g->nodes[i].dist_from_source < path_best.dist_cost) remaining++;
	  else g->nodes[i].active = false;
	  if (g->nodes[i].time_to_dest == 0) g->nodes[i].active = false;

	}

      //printf("Graph reduced: %u nodes remain active.\n", remaining);




      //printf("\nReducing graph forward using A* (distance) with budget W = %lf...\n", W);

      AStarResult astar_forward_reduction_distance = a_star_with_bound(g, s, -1, W, path_best.dist_cost, true, s, d, W, path_best.dist_cost, 1, -1); // forward, time

      remaining = 0;
      for (uint32_t i = 0; i < g->n_nodes; i++)
	if(g->nodes[i].active)
	if(i != s && i != d){
	  if (g->nodes[i].time_to_dest + g->nodes[i].time_from_source > 0 && g->nodes[i].time_to_dest + g->nodes[i].time_from_source <=  W && g->nodes[i].dist_to_dest + g->nodes[i].dist_from_source < path_best.dist_cost) remaining++;
	  else g->nodes[i].active = false;

	}

      //printf("Graph reduced: %u nodes remain active.\n", remaining);

      free_astar_result(&astar_forward_reduction_distance);

      //printf("\nReducing graph backward using A* (distance) with budget W = %lf...\n", W);

      AStarResult astar_backward_reduction_distance = a_star_with_bound(g, d, -1, W, path_best.dist_cost, false, s, d, W, path_best.dist_cost, 1, -1); // forward, time


      free_astar_result(&astar_backward_reduction_distance);


      // Conta i nodi rimasti

      remaining = 0;
      for (uint32_t i = 0; i < g->n_nodes; i++)
	if(g->nodes[i].active)
	if(i != s && i != d){
	  if (g->nodes[i].time_to_dest + g->nodes[i].time_from_source > 0 && g->nodes[i].time_to_dest + g->nodes[i].time_from_source <=  W && g->nodes[i].dist_to_dest + g->nodes[i].dist_from_source < path_best.dist_cost) remaining++;
	  else g->nodes[i].active = false;

	}



      //printf("Graph reduced: %u nodes remain active.\n", remaining);

      end_time_red = clock();

      //printf("\nStart INTERACTION heuristic...\n");


      clock_t start_exec_lambda = clock();
      clock_t end_exec_lambda = clock();
      clock_t start_time_best_lambda = clock();
      clock_t end_time_best_lambda = clock();
      double time_best_solution = 0.0;
      int nIt = 0;
      int max_iterations_red = INF;
      int iter_red = 0;
      int iter_lambda = 0;
      int best_iter_red = 0;
      int best_iter_lambda = 0;
      double best_lambda = 1;
      int check_connected = 1;
      while(path_inf.time_cost > W && ((double)((double)(end_exec_lambda - start_exec_lambda) / CLOCKS_PER_SEC)) < time_limit && iter_red < max_iterations_red){

	AStarResult astar_lambda = {0};

	if(iter_red > 0){
	  astar_lambda = a_star_with_bound(g, s, d, W, path_best.dist_cost, true, s, d, W, path_best.dist_cost, 0, -1); // forward, time

	  if(astar_lambda.time_cost == DBL_MAX){
	    //printf("Exit: Graph is no more connected.\n");
	    free_astar_result(&astar_lambda);

	    break;
	  }

	  if(astar_lambda.time_cost <= W)
	    if(astar_lambda.dist_cost < path_best.dist_cost){

	    end_exec_lambda = clock();
	    if(((double)((double)(end_exec_lambda - start_exec_lambda) / CLOCKS_PER_SEC)) > time_limit)
	      break;

	      free_astar_result(&path_best);
	      path_best = copy_res(&astar_lambda);
	      improvements_push(&improvements, path_best.dist_cost, path_best.time_cost,
				elapsed_seconds(start_exec), 0.0, nIt);

	      best_iter_red = iter_red;
	      best_iter_lambda = nIt;
	      best_lambda = 0;
	      end_time_best_lambda = clock();
	      end_exec = clock();
	      time_best_solution = ((double)((double)(end_exec - start_exec) / CLOCKS_PER_SEC));
	      /*
	      remaining = 0;
	      for (uint32_t i = 0; i < g->n_nodes; i++)
		if(i != s && i != d){
		  if (g->nodes[i].time_to_dest + g->nodes[i].time_from_source > 0 && g->nodes[i].time_to_dest + g->nodes[i].time_from_source <=  W && g->nodes[i].dist_to_dest + g->nodes[i].dist_from_source < path_best.dist_cost) remaining++;
		  else g->nodes[i].active = false;

		}
	      */
	    }

	  free_astar_result(&astar_lambda);
	}




	lambda = 0.5;

	double lambda_sup = 1, lambda_inf = 0;
	iter_lambda = 0;
	while(((double)((double)(end_exec_lambda - start_exec_lambda) / CLOCKS_PER_SEC)) < time_limit && iter_lambda < max_iterations){
	  //Call lambda A*
	  iter_lambda++;
	  nIt++;
	  free_astar_result(&astar_lambda);
	  clock_t start_exec_astar_lambda = clock();
	  astar_lambda = a_star_with_bound(g, s, d, W, path_best.dist_cost, true, s, d, W, path_best.dist_cost, lambda, -1); // forward, time
	  clock_t end_exec_astar_lambda = clock();

	  if(astar_lambda.time_cost <= W){

	    end_exec_lambda = clock();
	    if(((double)((double)(end_exec_lambda - start_exec_lambda) / CLOCKS_PER_SEC)) > time_limit){
	      check_connected = 0;

	      break;
	    }

	    if(astar_lambda.dist_cost < path_best.dist_cost){
	      free_astar_result(&path_best);
	      path_best = copy_res(&astar_lambda);
	      improvements_push(&improvements, path_best.dist_cost, path_best.time_cost,
				elapsed_seconds(start_exec), lambda, nIt);

	      best_iter_red = iter_red;
	      best_iter_lambda = nIt;
	      best_lambda = lambda;
	      end_time_best_lambda = clock();
	      end_exec = clock();
	      time_best_solution = ((double)((double)(end_exec - start_exec) / CLOCKS_PER_SEC));

	    }
	    lambda_inf = lambda;
	    lambda = lambda_sup - (lambda_sup - lambda_inf)/2;
	  }else{
	    lambda_sup = lambda;
	    lambda = lambda_sup - (lambda_sup - lambda_inf)/2;
	  }

	  end_exec_lambda = clock();
	}
	free_astar_result(&astar_lambda);


	if(check_connected == 0 || run_reduction_heuristic == 0)
	  break;


	free_astar_result(&path_inf);
	path_inf = a_star_with_bound(g, s, d, W, path_best.dist_cost, true, s, d, W, path_best.dist_cost, 1, -1); // forward, time


	if(path_inf.time_cost > W){
	  if(path_inf.time_cost == DBL_MAX){
	    //printf("Exit: Graph is no more connected.\n");

	    break;
	  }
	  W = path_inf.bound_path - 1e-06;

	  remaining = 0;
	  for (uint32_t i = 0; i < g->n_nodes; i++)
	    if(i != s && i != d){
	    if (g->nodes[i].active && g->nodes[i].time_to_dest + g->nodes[i].time_from_source > 0 && g->nodes[i].time_to_dest + g->nodes[i].time_from_source <  W && g->nodes[i].dist_to_dest + g->nodes[i].dist_from_source < path_best.dist_cost) remaining++;
	    else g->nodes[i].active = false;
	  }
	  //printf("Graph reduced: %u nodes remain active.\n", remaining);

	}else{
	  if(path_inf.dist_cost < path_best.dist_cost){

	    end_exec_lambda = clock();
	    if(((double)((double)(end_exec_lambda - start_exec_lambda) / CLOCKS_PER_SEC)) > time_limit)
	      break;

	    free_astar_result(&path_best);
	    path_best = copy_res(&path_inf);
	    improvements_push(&improvements, path_best.dist_cost, path_best.time_cost,
			      elapsed_seconds(start_exec), 1.0, nIt);

	    best_iter_red = iter_red;
	    best_iter_lambda = nIt;
	    best_lambda = 1;
	    end_time_best_lambda = clock();
	    end_exec = clock();
	    time_best_solution = ((double)((double)(end_exec - start_exec) / CLOCKS_PER_SEC));

	  }


	  break;
	}


	end_exec_lambda = clock();
	iter_red++;
      }
      end_exec_lambda = clock();
      end_exec = clock();

      //printf("Best solution found.\n");
      //printf("Cost (distance): %lf, Cost (time): %lf best_iter_red: %d best_iter_lambda %d best_lambda %lf time read graph: %.3f time-preprocessing: %.3f time heuristic %.3f total run time %.3f\n", path_best.dist_cost, path_best.time_cost, best_iter_red, best_iter_lambda, best_lambda, ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), ((double)((double)(end_time_red- start_time_red) / CLOCKS_PER_SEC)), ((double)((double)(end_exec_lambda- start_exec_lambda) / CLOCKS_PER_SEC)), ((double)((double)(end_exec- start_exec) / CLOCKS_PER_SEC)) + ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)));
      fprintf(fout, "%u %u %.0lf %.3f %.3f %.0f %.0f %.0f %.3f %.3f %.3f %.0f %.0f %d %lf %d %.3f %.3f\n", s, d, (double)W_read, ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), ((double)((double)(end_time_spt - start_time_spt) / CLOCKS_PER_SEC)), time_spt, dist_spt, dist_spd, ((double)((double)(end_time_red - start_time_red) / CLOCKS_PER_SEC)), ((double)((double)(end_exec_lambda- start_exec_lambda) / CLOCKS_PER_SEC)), ((double)((double)(end_time_best_lambda- start_time_best_lambda) / CLOCKS_PER_SEC)), path_best.time_cost, path_best.dist_cost, best_iter_lambda, best_lambda, nIt, ((double)((double)(end_exec- start_exec) / CLOCKS_PER_SEC)) + ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)), time_best_solution + ((double)((double)(end_read - start_read) / CLOCKS_PER_SEC)));

      fflush(fout);
      improvements_write_file(improvements_filename, &improvements);
      improvements_free(&improvements);
      // ---------- Cleanup ----------

      free_astar_result(&path_inf);
      free_astar_result(&path_best);



      if(runall != 1)
	return 0;


    }

    fclose(fout);
    freeGraph(g);


    if(runall == 1)
      fclose(input_runall);

    return 0;
}
