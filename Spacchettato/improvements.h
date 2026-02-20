#ifndef IMPROVEMENTS_H
#define IMPROVEMENTS_H

#include <stddef.h>
#include <time.h>

typedef struct {
    double length;
    double sunLength;
    double timeFound;
    double lambda;
    int iter;
} Improvement;

typedef struct {
    Improvement *items;
    size_t count;
    size_t capacity;
} ImprovementList;

void improvements_init(ImprovementList *list);
void improvements_free(ImprovementList *list);
void improvements_push(ImprovementList *list, double length, double sunLength,
                       double timeFound, double lambda, int iter);
void improvements_write_file(const char *filename, const ImprovementList *list);

double elapsed_seconds(clock_t start);

#endif
