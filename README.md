# matrix-multiplication-optimization

## Description
Optimizations for the matrix multiplication C = A * B.

## Compilation
gcc main.c mxm.c util.c -O3

## Implemented Methods
- Naive
- Matrix blocking
- Loop reordering
- Copying optimization
- Loop unrolling

## TODO
- Select an appropiate block size automatically
- Add SSE instructions
