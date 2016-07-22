# matrix-multiplication-optimization

## Description
Optimizations for the multiplication C = A * B of dense matrices.

## Compilation
gcc main.c mxm.c util.c -O3

## Implemented Methods
- Naive
- Matrix blocking
- Loop reordering
- Copying optimization
- Loop unrolling

## Speedup
- Machine: Intel Core i5 2500K
- Data: 2000x2000 matrices
- 7 against the -O3 naive version and 26 against the -O0 naive version

## TODO
- Select an appropiate block size automatically
- Add SSE instructions
