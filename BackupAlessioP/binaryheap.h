//
// Created by alessio on 26/05/25.
//

#ifndef BINARYHEAP_H
#define BINARYHEAP_H

// Heap node: stores the node index and the current distance (float).
typedef struct {
    int node;
    float dist;
} HeapNode;

// Min-heap structure.
typedef struct {
    int size;
    int capacity;
    int *pos;         // pos[v] stores the position of node v in the heap array.
    HeapNode* array;
} MinHeap;

// Creates a min-heap with the given capacity.
MinHeap* createMinHeap(int capacity) {
    MinHeap* heap = (MinHeap*)malloc(sizeof(MinHeap));
    heap->pos = (int*)malloc(capacity * sizeof(int));
    heap->size = 0;
    heap->capacity = capacity;
    heap->array = (HeapNode*)malloc(capacity * sizeof(HeapNode));
    return heap;
}

// Swaps two heap nodes.
void swapHeapNode(HeapNode* a, HeapNode* b) {
    HeapNode temp = *a;
    *a = *b;
    *b = temp;
}

// Restores the min-heap property starting from index idx.
void minHeapify(MinHeap* heap, int idx) {
    int smallest = idx;
    int left = 2 * idx + 1;
    int right = 2 * idx + 2;

    if (left < heap->size && heap->array[left].dist < heap->array[smallest].dist)
        smallest = left;
    if (right < heap->size && heap->array[right].dist < heap->array[smallest].dist)
        smallest = right;

    if (smallest != idx) {
        heap->pos[heap->array[smallest].node] = idx;
        heap->pos[heap->array[idx].node] = smallest;
        swapHeapNode(&heap->array[smallest], &heap->array[idx]);
        minHeapify(heap, smallest);
    }
}

// Checks if the heap is empty.
int isEmpty(MinHeap* heap) {
    return heap->size == 0;
}

// Extracts the node with the minimum distance from the heap.
int extractMin(MinHeap* heap) {
    if (isEmpty(heap))
        return -1;
    HeapNode root = heap->array[0];
    HeapNode lastNode = heap->array[heap->size - 1];
    heap->array[0] = lastNode;
    heap->pos[root.node] = heap->size - 1;
    heap->pos[lastNode.node] = 0;

    heap->size--;
    minHeapify(heap, 0);

    return root.node;
}

// Inserts a new node v with key f into the heap.
void insertHeap(MinHeap* h, int v, float f) {
    int i = h->size++;
    h->array[i].node = v;
    h->array[i].dist = f;
    h->pos[v] = i;
    // “Bubble up”
    while (i > 0 && h->array[i].dist < h->array[(i - 1) / 2].dist) {
        int p = (i - 1) / 2;
        h->pos[h->array[i].node] = p;
        h->pos[h->array[p].node] = i;
        swapHeapNode(&h->array[i], &h->array[p]);
        i = p;
    }
}

// Updates the distance of a node and reorders the heap.
void decreaseKey(MinHeap* heap, int node, float dist) {
    int i = heap->pos[node];
    heap->array[i].dist = dist;

    while (i && heap->array[i].dist < heap->array[(i - 1) / 2].dist) {
        heap->pos[heap->array[i].node] = (i - 1) / 2;
        heap->pos[heap->array[(i - 1) / 2].node] = i;
        swapHeapNode(&heap->array[i], &heap->array[(i - 1) / 2]);
        i = (i - 1) / 2;
    }
}

// Checks if a node is still present in the heap.
int isInMinHeap(MinHeap* heap, int node) {
    return (heap->pos[node] < heap->size);
}

// Frees the memory allocated for the heap.
void freeMinHeap(MinHeap* heap) {
    if (heap) {
        free(heap->pos);
        free(heap->array);
        free(heap);
    }
}

#endif // BINARYHEAP_H
