cmake_minimum_required(VERSION 3.28)
project(BinaryCombinedHeuristic C)

set(CMAKE_C_STANDARD 11)

file(COPY ${CMAKE_SOURCE_DIR}/data DESTINATION ${CMAKE_BINARY_DIR}) #Si copia all'interno della directory di build tutta la cartella data
# I risultati (gli output) si troveranno all'interno di /cmake-build-release(o debug)/data.

add_executable(BinaryCombinedHeuristic main.c
        graph.h
        Dijkstra.h
        BinaryCombinedHeuristic.h
        AStar.h
        spherical_heuristic.h
        binaryheap.h
)

target_link_libraries(BinaryCombinedHeuristic m)