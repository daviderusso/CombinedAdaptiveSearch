#include "improvements.h"

#include <stdio.h>
#include <stdlib.h>

void improvements_init(ImprovementList *list) {
    list->items = NULL;
    list->count = 0;
    list->capacity = 0;
}

void improvements_free(ImprovementList *list) {
    free(list->items);
    list->items = NULL;
    list->count = 0;
    list->capacity = 0;
}

void improvements_push(ImprovementList *list, double length, double sunLength,
                       double timeFound, double lambda, int iter) {
    if (list->count == list->capacity) {
        size_t new_cap = list->capacity ? list->capacity * 2 : 8;
        Improvement *tmp = realloc(list->items, new_cap * sizeof(Improvement));
        if (!tmp) exit(EXIT_FAILURE);
        list->items = tmp;
        list->capacity = new_cap;
    }
    list->items[list->count++] = (Improvement){length, sunLength, timeFound, lambda, iter};
}

void improvements_write_file(const char *filename, const ImprovementList *list) {
    FILE *f = fopen(filename, "w");
    if (!f) {
        perror("Error opening improvements file");
        return;
    }
    fprintf(f, "length sunLength timeFound lambda iter\n");
    for (size_t i = 0; i < list->count; i++) {
        const Improvement *imp = &list->items[i];
        fprintf(f, "%.0f %.0f %.6f %.6f %d\n",
                imp->length, imp->sunLength, imp->timeFound, imp->lambda, imp->iter);
    }
    fclose(f);
}

double elapsed_seconds(clock_t start) {
    return ((double) (clock() - start)) / CLOCKS_PER_SEC;
}
