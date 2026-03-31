#include <stdio.h> // for printf
#include <stdlib.h> // for rand and srand
#include <time.h> // for time

#define SIZE 5 // define the size of the square matrices

// in this task, we compose a dense square matrix multiplication program in C. Generate two random 5x5 matrices and compute six matrix products using six different (i,j,k) orders. 



/* 
visual representation of matrix multiplication:
------------------------------     ----------------------------      -----------------------------
|   a1    a2    a3    a4    a5 |   | b1    b2    b3    b4    b5 |     | c1    c2    c3    c4    c5 |
|   a6    a7    a8    a9    a10|   | b6    b7    b8    b9    b10|     | c6    c7    c8    c9    c10|
|   a11   a12   a13   a14   a15| * | b11   b12   b13   b14   b15|  =  | c11   c12   c13   c14   c15|
|   a16   a17   a18   a19   a20|   | b16   b17   b18   b19   b20|     | c16   c17   c18   c19   c20|
|   a21   a22   a23   a24   a25|   | b21   b22   b23   b24   b25|     | c21   c22   c23   c24   c25|
----------------------------        -----------------------------     -----------------------------
*/

// matrix i j k

void matrix_multiply_ijk(double A[SIZE][SIZE], double B[SIZE][SIZE], double C[SIZE][SIZE]){ // A, B, C are 5x5 matrices
    for (int i = 0; i < SIZE; i++) // loop over rows of A and C
        for (int j = 0; j < SIZE; j++){ // loop over columns of B and C
            double sum = 0.0; // initialize sum for C[i][j]. we could also initialize C[i][j] to 0.0 before the innermost loop, but using a separate sum variable allows us to accumulate the product of A's row and B's column without modifying C until we have the final sum, which can be more efficient.
            for (int k = 0; k < SIZE; k++) // loop over columns of A and rows of B
                sum += A[i][k] * B[k][j]; // accumulate the product of A's row and B's column into sum
            C[i][j] += sum; // add the computed sum to C[i][j] 
        }   
}

// this gives us something like this:
// a1 * b1 + a2 * b6 + a3 * b11 + a4 * b16 + a5 * b21 = c1      --> i = 0, j = 0, k = 0,1,2,3,4
// a1 * b2 + a2 * b7 + a3 * b12 + a4 * b17 + a5 * b22 = c2      --> i = 0, j = 1, k = 0,1,2,3,4
// a1 * b3 + a2 * b8 + a3 * b13 + a4 * b18 + a5 * b23 = c3      --> i = 0, j = 2, k = 0,1,2,3,4
// a1 * b4 + a2 * b9 + a3 * b14 + a4 * b19 + a5 * b24 = c4      --> i = 0, j = 3, k = 0,1,2,3,4
// a1 * b5 + a2 * b10 + a3 * b15 + a4 * b20 + a5 * b25 = c5     --> i = 0, j = 4, k = 0,1,2,3,4



// matrix j i k 

void matrix_multiply_jik(double A[SIZE][SIZE], double B[SIZE][SIZE], double C[SIZE][SIZE]){ // A, B, C are 5x5 matrices
    for (int j = 0; j < SIZE; j++) // loop over columns of B and C
        for (int i = 0; i < SIZE; i++){ // loop over rows of A and C
            double sum = 0.0; // initialize sum for C[i][j]. we could also initialize C[i][j] to 0.0 before the innermost loop, but using a separate sum variable allows us to accumulate the product of A's row and B's column without modifying C until we have the final sum, which can be more efficient.
            for (int k = 0; k < SIZE; k++) // loop over columns of A and rows of B
                sum += A[i][k] * B[k][j]; // accumulate the product of A's row and B's column into sum
            C[i][j] += sum; // add the computed sum to C[i][j]
        }
}

// this gives us something like this:
// a1 * b1 + a2 * b6 + a3 * b11 + a4 * b16 + a5 * b21 = c1          --> j = 0, i = 0, k = 0,1,2,3,4
// a6 * b1 + a7 * b6 + a8 * b11 + a9 * b16 + a10 * b21 = c6         --> j = 0, i = 1, k = 0,1,2,3,4
// a11 * b1 + a12 * b6 + a13 * b11 + a14 * b16 + a15 * b21 = c11    --> j = 0, i = 2, k = 0,1,2,3,4
// a16 * b1 + a17 * b6 + a18 * b11 + a19 * b16 + a20 * b21 = c16    --> j = 0, i = 3, k = 0,1,2,3,4
// a21 * b1 + a22 * b6 + a23 * b11 + a24 * b16 + a25 * b21 = c21    --> j = 0, i = 4, k = 0,1,2,3,4


// matrix j k i

void matrix_multiply_jki(double A[SIZE][SIZE], double B[SIZE][SIZE], double C[SIZE][SIZE]){ // A, B, C are 5x5 matrices
    double r; // temporary variable to hold A[k][j]. this helps to reduce the number of memory accesses to A, which can improve cache performance.
    for (int j = 0; j < SIZE; j++) // loop over columns of B and C
        for (int k = 0; k < SIZE; k++){ // loop over columns of A and rows of B
            r = B[k][j]; // store B[k][j] in r to avoid multiple memory accesses to B[k][j] in the innermost loop
            for (int i = 0; i < SIZE; i++) // loop over rows of A and C
                C[i][j] += A[i][k] * r; // accumulate the product of A's row and B's column into C[i][j]. using r instead of B[k][j] can improve cache performance by reducing memory accesses to B.
        }
}

// this gives us something like this:
// a1 * b1 + a6 * b6 + a11 * b11 + a16 * b16 + a21 * b21 = c1          --> j = 0, k = 0, i = 0,1,2,3,4
// a1 * b2 + a6 * b7 + a11 * b12 + a16 * b17 + a21 * b22 = c2          --> j = 1, k = 0, i = 0,1,2,3,4
// a1 * b3 + a6 * b8 + a11 * b13 + a16 * b18 + a21 * b23 = c3          --> j = 2, k = 0, i = 0,1,2,3,4
// a1 * b4 + a6 * b9 + a11 * b14 + a16 * b19 + a21 * b24 = c4          --> j = 3, k = 0, i = 0,1,2,3,4
// a1 * b5 + a6 * b10 + a11 * b15 + a16 * b20 + a21 * b25 = c5         --> j = 4, k = 0, i = 0,1,2,3,4

void matrix_multiply_kji(double A[SIZE][SIZE], double B[SIZE][SIZE], double C[SIZE][SIZE]){ // A, B, C are 5x5 matrices
    double r; // temporary variable to hold A[k][j]. this helps to reduce the number of memory accesses to A, which can improve cache performance.
    for (int k = 0; k < SIZE; k++) // loop over columns of A and rows of B
        for (int j = 0; j < SIZE; j++){ // loop over columns of B and C
            r = B[k][j]; // store B[k][j] in r to avoid multiple memory accesses to B[k][j] in the innermost loop
            for (int i = 0; i < SIZE; i++) // loop over rows of A and C
                C[i][j] += A[i][k] * r; // accumulate the product of A's row and B's column into C[i][j]. using r instead of B[k][j] can improve cache performance by reducing memory accesses to B.
        }
}


// this gives us something like this:
// a1 * b1 + a6 * b6 + a11 * b11 + a16 * b16 + a21 * b21 = c1          --> k = 0, j = 0, i = 0,1,2,3,4
// a1 * b2 + a6 * b7 + a11 * b12 + a16 * b17 + a21 * b22 = c2          --> k = 0, j = 1, i = 0,1,2,3,4
// a1 * b3 + a6 * b8 + a11 * b13 + a16 * b18 + a21 * b23 = c3          --> k = 0, j = 2, i = 0,1,2,3,4
// a1 * b4 + a6 * b9 + a11 * b14 + a16 * b19 + a21 * b24 = c4          --> k = 0, j = 3, i = 0,1,2,3,4
// a1 * b5 + a6 * b10 + a11 * b15 + a16 * b20 + a21 * b25 = c5         --> k = 0, j = 4, i = 0,1,2,3,4





void matrix_multiply_kij(double A[SIZE][SIZE], double B[SIZE][SIZE], double C[SIZE][SIZE]){ // A, B, C are 5x5 matrices
    double r;  // temporary variable to hold A[k][j]. this helps to reduce the number of memory accesses to A, which can improve cache performance.
    for (int k = 0; k < SIZE; k++) // loop over columns of A and rows of B
        for (int i = 0; i < SIZE; i++){ // loop over rows of A and C
            r = A[i][k]; // store A[i][k] in r to avoid multiple memory accesses to A[i][k] in the innermost loop
            for (int j = 0; j < SIZE; j++) // loop over columns of B and C
                C[i][j] += r * B[k][j]; // accumulate the product of A's row and B's column into C[i][j]. using r instead of A[i][k] can improve cache performance by reducing memory accesses to A.
        }
}

// this gives us something like this:
// a1 * b1 + a6 * b6 + a11 * b11 + a16 * b16 + a21 * b21 = c1          --> k = 0, i = 0, j = 0,1,2,3,4
// a1 * b2 + a6 * b7 + a11 * b12 + a16 * b17 + a21 * b22 = c2          --> k = 0, i = 0, j = 1,1,2,3,4
// a1 * b3 + a6 * b8 + a11 * b13 + a16 * b18 + a21 * b23 = c3          --> k = 0, i = 0, j = 2,1,2,3,4
// a1 * b4 + a6 * b9 + a11 * b14 + a16 * b19 + a21 * b24 = c4          --> k = 0, i = 0, j = 3,1,2,3,4
// a1 * b5 + a6 * b10 + a11 * b15 + a16 * b20 + a21 * b25 = c5         --> k = 0, i = 0, j = 4,1,2,3,4


void matrix_multiply_ikj(double A[SIZE][SIZE], double B[SIZE][SIZE], double C[SIZE][SIZE]){ // A, B, C are 5x5 matrices
    double r; // temporary variable to hold A[i][k]. this helps to reduce the number of memory accesses to A, which can improve cache performance.
    for (int i = 0; i < SIZE; i++) // loop over rows of A and C
        for (int k = 0; k < SIZE; k++){ // loop over columns of A and rows of B
            r = A[i][k]; // store A[i][k] in r to avoid multiple memory accesses to A[i][k] in the innermost loop
            for (int j = 0; j < SIZE; j++) // loop over columns of B and C
                C[i][j] += r * B[k][j]; // accumulate the product of A's row and B's column into C[i][j]. using r instead of A[i][k] can improve cache performance by reducing memory accesses to A.
        }
}


// this gives us something like this:
// a1 * b1 + a6 * b6 + a11 * b11 + a16 * b16 + a21 * b21 = c1          --> i = 0, k = 0, j = 0,1,2,3,4
// a1 * b2 + a6 * b7 + a11 * b12 + a16 * b17 + a21 * b22 = c2          --> i = 0, k = 0, j = 1,1,2,3,4
// a1 * b3 + a6 * b8 + a11 * b13 + a16 * b18 + a21 * b23 = c3          --> i = 0, k = 0, j = 2,1,2,3,4
// a1 * b4 + a6 * b9 + a11 * b14 + a16 * b19 + a21 * b24 = c4          --> i = 0, k = 0, j = 3,1,2,3,4
// a1 * b5 + a6 * b10 + a11 * b15 + a16 * b20 + a21 * b25 = c5         --> i = 0, k = 0, j = 4,1,2,3,4

// Function to print a matrix
void print_matrix(double C[SIZE][SIZE]) // this function takes a 5x5 matrix C as input and prints its elements in a formatted way. it uses nested loops to iterate through the rows and columns of the matrix, printing each element with a specific width for alignment. after printing each row, it prints a newline character to move to the next line.
{
    for (int i = 0; i < SIZE; i++) // loop over rows of C
    {
        for (int j = 0; j < SIZE; j++) // loop over columns of C
        {
            printf("%5.1f ", C[i][j]); // print each element of C with a width of 5 characters and 1 decimal place, followed by a space. this formatting ensures that the elements are aligned in columns when printed.
        }
        printf("\n"); // print a newline character after each row to move to the next line for the next row of the matrix.
    }
}





int main()
{
    // Seed random number generator
    srand(time(NULL)); // seeding the random number generator with the current time ensures that we get different random values each time we run the program, which is important for testing the matrix multiplication functions with different inputs.

    double A[SIZE][SIZE], B[SIZE][SIZE]; // Create 2 input matrices 

    // Create 6 result matrices 
    double C1[SIZE][SIZE];
    double C2[SIZE][SIZE];
    double C3[SIZE][SIZE];
    double C4[SIZE][SIZE];
    double C5[SIZE][SIZE];
    double C6[SIZE][SIZE];

    // Initialize matrices A and B with random values
    for (int i = 0; i < SIZE; i++) // loop over rows of A and C
    {
        for (int j = 0; j < SIZE; j++) // loop over columns of B and C
        {
            A[i][j] = rand() % 10; // initialize A[i][j] with a random integer between 0 and 9
            B[i][j] = rand() % 10; // initialize B[i][j] with a random integer between 0 and 9

            // Initialize all result matrices to 0
            C1[i][j] = 0.0;
            C2[i][j] = 0.0;
            C3[i][j] = 0.0;
            C4[i][j] = 0.0;
            C5[i][j] = 0.0;
            C6[i][j] = 0.0;
        }
    }

    // Call all 6 multiplication functions
    matrix_multiply_ijk(A, B, C1); 
    matrix_multiply_jik(A, B, C2);
    matrix_multiply_jki(A, B, C3);
    matrix_multiply_kji(A, B, C4);
    matrix_multiply_kij(A, B, C5);
    matrix_multiply_ikj(A, B, C6);



    // Print all result matrices
    printf("\nResult (ijk):\n");
    print_matrix(C1);

    printf("\nResult (jik):\n");
    print_matrix(C2);

    printf("\nResult (jki):\n");
    print_matrix(C3);

    printf("\nResult (kji):\n");
    print_matrix(C4);

    printf("\nResult (kij):\n");
    print_matrix(C5);

    printf("\nResult (ikj):\n");
    print_matrix(C6);

    return 0;
}
