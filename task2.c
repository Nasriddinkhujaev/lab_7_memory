#include <stdio.h>
#include <stdlib.h>
#include <time.h>

// ---------- MATRIX MULTIPLICATION FUNCTIONS ----------

void matrix_multiply_ijk(int n, double A[n][n], double B[n][n], double C[n][n]) // A, B, C are n x n matrices
{
    for (int i = 0; i < n; i++) // loop over rows of A and C
        for (int j = 0; j < n; j++) // loop over columns of B and C
        {
            double sum = 0.0; // temporary variable to hold the sum of products for C[i][j]. this helps to reduce the number of memory accesses to C[i][j], which can improve cache performance.
            for (int k = 0; k < n; k++) // loop over columns of A and rows of B
                sum += A[i][k] * B[k][j]; // accumulate the product of A's row and B's column into sum. this is the standard way to compute the dot product of A's row and B's column for C[i][j].
            C[i][j] = sum; // store the computed value of C[i][j] after the innermost loop. this ensures that we only write to C[i][j] once per element, which can improve cache performance by reducing memory writes.
        }
}

void matrix_multiply_jik(int n, double A[n][n], double B[n][n], double C[n][n]) // A, B, C are n x n matrices
{
    for (int j = 0; j < n; j++) // loop over columns of B and C
        for (int i = 0; i < n; i++) // loop over rows of A and C
        {
            double sum = 0.0; // temporary variable to hold the sum of products for C[i][j]. this helps to reduce the number of memory accesses to C[i][j], which can improve cache performance.
            for (int k = 0; k < n; k++) // loop over columns of A and rows of B
                sum += A[i][k] * B[k][j]; // accumulate the product of A's row and B's column into sum. this is the standard way to compute the dot product of A's row and B's column for C[i][j].
            C[i][j] = sum; // store the computed value of C[i][j] after the innermost loop. this ensures that we only write to C[i][j] once per element, which can improve cache performance by reducing memory writes.
        }
}

void matrix_multiply_jki(int n, double A[n][n], double B[n][n], double C[n][n]) // A, B, C are n x n matrices
{
    for (int j = 0; j < n; j++) // loop over columns of B and C
        for (int k = 0; k < n; k++) // loop over columns of A and rows of B
        {
            double r = B[k][j]; // store B[k][j] in r to avoid multiple memory accesses to B[k][j] in the innermost loop. this can improve cache performance by reducing memory accesses to B.
            for (int i = 0; i < n; i++) // loop over rows of A and C
                C[i][j] += A[i][k] * r; // accumulate the product of A's row and B's column into C[i][j]. using r instead of B[k][j] can improve cache performance by reducing memory accesses to B.
        }
}

void matrix_multiply_kji(int n, double A[n][n], double B[n][n], double C[n][n]) // A, B, C are n x n matrices
{
    for (int k = 0; k < n; k++) // loop over columns of A and rows of B
        for (int j = 0; j < n; j++) // loop over columns of B and C
        {
            double r = B[k][j]; // store B[k][j] in r to avoid multiple memory accesses to B[k][j] in the innermost loop. this can improve cache performance by reducing memory accesses to B.
            for (int i = 0; i < n; i++) // loop over rows of A and C
                C[i][j] += A[i][k] * r; // accumulate the product of A's row and B's column into C[i][j]. using r instead of B[k][j] can improve cache performance by reducing memory accesses to B.
        }
}

void matrix_multiply_kij(int n, double A[n][n], double B[n][n], double C[n][n]) // A, B, C are n x n matrices
{
    for (int k = 0; k < n; k++) // loop over columns of A and rows of B
        for (int i = 0; i < n; i++) // loop over rows of A and C
        {
            double r = A[i][k]; // store A[i][k] in r to avoid multiple memory accesses to A[i][k] in the innermost loop. this can improve cache performance by reducing memory accesses to A.
            for (int j = 0; j < n; j++) // loop over columns of B and C
                C[i][j] += r * B[k][j]; // accumulate the product of A's row and B's column into C[i][j]. using r instead of A[i][k] can improve cache performance by reducing memory accesses to A.
        }
}

void matrix_multiply_ikj(int n, double A[n][n], double B[n][n], double C[n][n]) // A, B, C are n x n matrices
{
    for (int i = 0; i < n; i++) // loop over rows of A and C
        for (int k = 0; k < n; k++) // loop over columns of A and rows of B
        {
            double r = A[i][k]; // store A[i][k] in r to avoid multiple memory accesses to A[i][k] in the innermost loop. this can improve cache performance by reducing memory accesses to A.
            for (int j = 0; j < n; j++) // loop over columns of B and C
                C[i][j] += r * B[k][j]; // accumulate the product of A's row and B's column into C[i][j]. using r instead of A[i][k] can improve cache performance by reducing memory accesses to A.
        }
}


// ---------- MAIN ----------


int main()
{

    srand(time(NULL)); // seed the random number generator with the current time to ensure different random values each time the program is run.

    int sizes[] = {20, 40, 80, 160, 320, 640, 1280}; // array of matrix sizes to test. 
    int num_sizes = 7; // number of different matrix sizes to test.

    printf("Size\tijk\tjik\tjki\tkji\tkij\tikj\n"); // print the header for the output table. this will help us to understand which column corresponds to which loop order when we print the results of the matrix multiplication timings.

    for (int s = 0; s < num_sizes; s++) // loop over different matrix sizes to test. for each size, we will allocate matrices, initialize them, perform the multiplications, and print the results.
    {

        int n = sizes[s]; // current matrix size to test. this will be used to allocate n x n matrices and to call the multiplication functions with the correct size.

        // Allocate matrices
        double (*A)[n] = malloc(sizeof(double[n][n]));
        double (*B)[n] = malloc(sizeof(double[n][n]));
        double (*C1)[n] = malloc(sizeof(double[n][n]));
        double (*C2)[n] = malloc(sizeof(double[n][n]));
        double (*C3)[n] = malloc(sizeof(double[n][n]));
        double (*C4)[n] = malloc(sizeof(double[n][n]));
        double (*C5)[n] = malloc(sizeof(double[n][n]));
        double (*C6)[n] = malloc(sizeof(double[n][n]));

        // Initialize matrices
        for (int i = 0; i < n; i++) // loop over rows of A and C
            for (int j = 0; j < n; j++) // loop over columns of B and C
            {
                A[i][j] = rand() % 10; // initialize A[i][j] with a random integer between 0 and 9
                B[i][j] = rand() % 10; // initialize B[i][j] with a random integer between 0 and 9

                C1[i][j] = C2[i][j] = C3[i][j] = 0.0; // initialize all result matrices to 0. this is important because the multiplication functions accumulate values into C, so we need to start with 0 for each element.
                C4[i][j] = C5[i][j] = C6[i][j] = 0.0; // initialize all result matrices to 0. this is important because the multiplication functions accumulate values into C, so we need to start with 0 for each element.
            }

        clock_t start, end; // variables to hold the start and end times for measuring the execution time of each multiplication function. we will use these to calculate the time taken for each loop order and print the results in a formatted way.

        // ijk
        start = clock(); // record the start time before calling the matrix multiplication function for the ijk loop order.
        matrix_multiply_ijk(n, A, B, C1); // call the matrix multiplication function for the ijk loop order with the current matrix size n and the allocated matrices A, B, and C1 to store the result.
        end = clock(); // record the end time after the function call to calculate the time taken for the ijk loop order.
        double t1 = (double)(end - start) / CLOCKS_PER_SEC; // calculate the time taken for the ijk loop order by taking the difference between end and start times and dividing by CLOCKS_PER_SEC to convert to seconds.

        // jik
        start = clock();
        matrix_multiply_jik(n, A, B, C2);
        end = clock();
        double t2 = (double)(end - start) / CLOCKS_PER_SEC;

        // jki
        start = clock();
        matrix_multiply_jki(n, A, B, C3);
        end = clock();
        double t3 = (double)(end - start) / CLOCKS_PER_SEC;

        // kji
        start = clock();
        matrix_multiply_kji(n, A, B, C4);
        end = clock();
        double t4 = (double)(end - start) / CLOCKS_PER_SEC;

        // kij
        start = clock();
        matrix_multiply_kij(n, A, B, C5);
        end = clock();
        double t5 = (double)(end - start) / CLOCKS_PER_SEC;

        // ikj
        start = clock();
        matrix_multiply_ikj(n, A, B, C6);
        end = clock();
        double t6 = (double)(end - start) / CLOCKS_PER_SEC;

        // Print results
        printf("%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\n",  
                n, t1, t2, t3, t4, t5, t6); // print the results for the current matrix size n and the time taken for each loop order in a formatted way. this will allow us to compare the performance of different loop orders for matrix multiplication as the size of the matrices increases.

        free(A); // free the allocated memory for matrix A to avoid memory leaks. this is important because we are allocating new matrices for each size in the loop, so we need to free the previous ones before moving on to the next size.
        free(B);
        free(C1); // free the allocated memory for matrix C1 to avoid memory leaks. this is important because we are allocating new matrices for each size in the loop, so we need to free the previous ones before moving on to the next size.
        free(C2);
        free(C3);
        free(C4);
        free(C5);
        free(C6);
    }

    return 0;
}