#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// #define SIZE1 20;
// #define SIZE2 40;
// #define SIZE3 80;
// #define SIZE4 160;
// #define SIZE5 320;
// #define SIZE6 640;
// #define SIZE7 1280;

// in this task, we generate random dense square matrices, with size of [20, 40, 80, 160, 320, 640, 1280]. Record and plot elapsed time using six different (i, j, k) orders. 
int size;

void matrix_multiply_ijk(double A[size][size], double B[size][size], double C[size][size]){ // A, B, C are 5x5 matrices
    for (int i = 0; i < size; i++) // loop over rows of A and C
        for (int j = 0; j < size; j++){ // loop over columns of B and C
            double sum = 0.0; // initialize sum for C[i][j]. we could also initialize C[i][j] to 0.0 before the innermost loop, but using a separate sum variable allows us to accumulate the product of A's row and B's column without modifying C until we have the final sum, which can be more efficient.
            for (int k = 0; k < size; k++) // loop over columns of A and rows of B
                sum += A[i][k] * B[k][j]; // accumulate the product of A's row and B's column into sum
            C[i][j] += sum; // add the computed sum to C[i][j] 
        }   
}

void matrix_multiply_jik(double A[size][size], double B[size][size], double C[size][size]){ // A, B, C are 5x5 matrices
    for (int j = 0; j < size; j++) // loop over columns of B and C
        for (int i = 0; i < size; i++){ // loop over rows of A and C
            double sum = 0.0; // initialize sum for C[i][j]. we could also initialize C[i][j] to 0.0 before the innermost loop, but using a separate sum variable allows us to accumulate the product of A's row and B's column without modifying C until we have the final sum, which can be more efficient.
            for (int k = 0; k < size; k++) // loop over columns of A and rows of B
                sum += A[i][k] * B[k][j]; // accumulate the product of A's row and B's column into sum
            C[i][j] += sum; // add the computed sum to C[i][j]
        }
}

void matrix_multiply_jki(double A[size][size], double B[size][size], double C[size][size]){ // A, B, C are 5x5 matrices
    double r; // temporary variable to hold A[k][j]. this helps to reduce the number of memory accesses to A, which can improve cache performance.
    for (int j = 0; j < size; j++) // loop over columns of B and C
        for (int k = 0; k < size; k++){ // loop over columns of A and rows of B
            r = B[k][j]; // store B[k][j] in r to avoid multiple memory accesses to B[k][j] in the innermost loop
            for (int i = 0; i < size; i++) // loop over rows of A and C
                C[i][j] += A[i][k] * r; // accumulate the product of A's row and B's column into C[i][j]. using r instead of B[k][j] can improve cache performance by reducing memory accesses to B.
        }
}

void matrix_multiply_kji(double A[size][size], double B[size][size], double C[size][size]){ // A, B, C are 5x5 matrices
    double r; // temporary variable to hold A[k][j]. this helps to reduce the number of memory accesses to A, which can improve cache performance.
    for (int k = 0; k < size; k++) // loop over columns of A and rows of B
        for (int j = 0; j < size; j++){ // loop over columns of B and C
            r = B[k][j]; // store B[k][j] in r to avoid multiple memory accesses to B[k][j] in the innermost loop
            for (int i = 0; i < size; i++) // loop over rows of A and C
                C[i][j] += A[i][k] * r; // accumulate the product of A's row and B's column into C[i][j]. using r instead of B[k][j] can improve cache performance by reducing memory accesses to B.
        }
}

void matrix_multiply_kij(double A[size][size], double B[size][size], double C[size][size]){ // A, B, C are 5x5 matrices
    double r;  // temporary variable to hold A[k][j]. this helps to reduce the number of memory accesses to A, which can improve cache performance.
    for (int k = 0; k < size; k++) // loop over columns of A and rows of B
        for (int i = 0; i < size; i++){ // loop over rows of A and C
            r = A[i][k]; // store A[i][k] in r to avoid multiple memory accesses to A[i][k] in the innermost loop
            for (int j = 0; j < size; j++) // loop over columns of B and C
                C[i][j] += r * B[k][j]; // accumulate the product of A's row and B's column into C[i][j]. using r instead of A[i][k] can improve cache performance by reducing memory accesses to A.
        }
}   

void matrix_multiply_ikj(double A[size][size], double B[size][size], double C[size][size]){ // A, B, C are 5x5 matrices
    double r; // temporary variable to hold A[i][k]. this helps to reduce the number of memory accesses to A, which can improve cache performance.
    for (int i = 0; i < size; i++) // loop over rows of A and C
        for (int k = 0; k < size; k++){ // loop over columns of A and rows of B
            r = A[i][k]; // store A[i][k] in r to avoid multiple memory accesses to A[i][k] in the innermost loop
            for (int j = 0; j < size; j++) // loop over columns of B and C
                C[i][j] += r * B[k][j]; // accumulate the product of A's row and B's column into C[i][j]. using r instead of A[i][k] can improve cache performance by reducing memory accesses to A.
        }
}

int main(){
    
    return 0;
}
