#ifndef HEAP_H
#define HEAP_H

#include <stdint.h>
#include <stdbool.h>

typedef struct {
    uint32_t node;
    double cost;
} HeapNode;

typedef struct {
    HeapNode *data;
    uint32_t size;
    uint32_t capacity;
} MinHeap;

MinHeap *createHeap(uint32_t capacity);
void push(MinHeap *h, uint32_t node, double cost);
HeapNode pop(MinHeap *h);
bool isEmpty(MinHeap *h);

#endif
