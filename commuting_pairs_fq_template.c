#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>

#define N 2
#define MAX_Q 11  // Allow q up to 11
int Q = 2;

int add(int a, int b) { return (a + b) % Q; }
int mul(int a, int b) { return (a * b) % Q; }

void matmul(int A[4], int B[4], int C[4]) {
    C[0] = add(mul(A[0], B[0]), mul(A[1], B[2]));
    C[1] = add(mul(A[0], B[1]), mul(A[1], B[3]));
    C[2] = add(mul(A[2], B[0]), mul(A[3], B[2]));
    C[3] = add(mul(A[2], B[1]), mul(A[3], B[3]));
}

bool commute(int A[4], int B[4]) {
    int AB[4], BA[4];
    matmul(A, B, AB);
    matmul(B, A, BA);
    for (int i = 0; i < 4; i++) if (AB[i] != BA[i]) return false;
    return true;
}

void int_to_matrix(int idx, int M[4]) {
    for (int i = 0; i < 4; i++) {
        M[i] = idx % Q;
        idx /= Q;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        printf("Usage: %s <prime_q>\n", argv[0]);
        return 1;
    }
    Q = atoi(argv[1]);
    int total = 0;
    int max_idx = (int)(pow(Q, 4));
    int A[4], B[4];

    for (int i = 0; i < max_idx; i++) {
        int_to_matrix(i, A);
        for (int j = 0; j < max_idx; j++) {
            int_to_matrix(j, B);
            if (commute(A, B)) total++;
        }
    }

    printf("%d\n", total);
    return 0;
}
