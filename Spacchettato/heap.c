#include "heap.h"

#include <stdlib.h>
#include <float.h>

static void swap(HeapNode *a, HeapNode *b) {
    HeapNode t = *a;
    *a = *b;
    *b = t;
}

static void heapifyUp(MinHeap *h, uint32_t i) {
    while (i && h->data[i].cost < h->data[(i - 1) / 2].cost) {
        swap(&h->data[i], &h->data[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

static void heapifyDown(MinHeap *h, uint32_t i) {
    uint32_t smallest = i;
    uint32_t l = 2 * i + 1, r = 2 * i + 2;
    if (l < h->size && h->data[l].cost < h->data[smallest].cost) smallest = l;
    if (r < h->size && h->data[r].cost < h->data[smallest].cost) smallest = r;
    if (smallest != i) {
        swap(&h->data[i], &h->data[smallest]);
        heapifyDown(h, smallest);
    }
}

MinHeap *createHeap(uint32_t capacity) {
    MinHeap *h = malloc(sizeof(MinHeap));
    if (!h) exit(EXIT_FAILURE);
    h->data = malloc(sizeof(HeapNode) * (capacity ? capacity : 4));
    if (!h->data) exit(EXIT_FAILURE);
    h->size = 0;
    h->capacity = capacity ? capacity : 4;
    return h;
}

void push(MinHeap *h, uint32_t node, double cost) {
    if (h->size == h->capacity) {
        h->capacity = h->capacity ? h->capacity * 2 : 4;
        HeapNode *tmp = realloc(h->data, sizeof(HeapNode) * h->capacity);
        if (!tmp) exit(EXIT_FAILURE);
        h->data = tmp;
    }
    h->data[h->size++] = (HeapNode){node, cost};
    heapifyUp(h, h->size - 1);
}

HeapNode pop(MinHeap *h) {
    if (h->size == 0) return (HeapNode){UINT32_MAX, DBL_MAX};
    HeapNode root = h->data[0];
    h->data[0] = h->data[--h->size];
    heapifyDown(h, 0);
    return root;
}

bool isEmpty(MinHeap *h) { return h->size == 0; }
