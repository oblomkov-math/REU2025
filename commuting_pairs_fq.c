#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define N 2         // Matrix size
#define Q 2        // Field size (set to 2, 3, 4, or 5)
#define FIELD_SIZE Q
#define MAT_ENTRIES (N * N)
#define NUM_MATRICES (int)(pow(Q, MAT_ENTRIES))

// Modular addition and multiplication
int add(int a, int b) {
    return (a + b) % Q;
}

int mul(int a, int b) {
    return (a * b) % Q;
}

// Multiply 2x2 matrices A and B, store result in C
void matmul(int A[MAT_ENTRIES], int B[MAT_ENTRIES], int C[MAT_ENTRIES]) {
    C[0] = add(mul(A[0], B[0]), mul(A[1], B[2]));
    C[1] = add(mul(A[0], B[1]), mul(A[1], B[3]));
    C[2] = add(mul(A[2], B[0]), mul(A[3], B[2]));
    C[3] = add(mul(A[2], B[1]), mul(A[3], B[3]));
}

// Check if A and B commute
bool commute(int A[MAT_ENTRIES], int B[MAT_ENTRIES]) {
    int AB[MAT_ENTRIES], BA[MAT_ENTRIES];
    matmul(A, B, AB);
    matmul(B, A, BA);
    for (int i = 0; i < MAT_ENTRIES; i++) {
        if (AB[i] != BA[i]) return false;
    }
    return true;
}

// Decode integer to matrix over F_q
void int_to_matrix(int idx, int M[MAT_ENTRIES]) {
    for (int i = 0; i < MAT_ENTRIES; i++) {
        M[i] = idx % Q;
        idx /= Q;
    }
}

int main() {
    int total = 0;
    int A[MAT_ENTRIES], B[MAT_ENTRIES];

    for (int i = 0; i < NUM_MATRICES; i++) {
        int_to_matrix(i, A);
        for (int j = 0; j < NUM_MATRICES; j++) {
            int_to_matrix(j, B);
            if (commute(A, B)) {
                total++;
            }
        }
    }

    printf("Total commuting pairs for 2x2 matrices over F_%d: %d\n", Q, total);
    return 0;
}


