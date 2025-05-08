#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define N 10
#define MAX_ITER 100
#define EPSILON 1e-6

void print_matrix(double A[N][N], int n, FILE *f) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%10.6lf ", A[i][j]);
            if (f) fprintf(f, "%10.6lf ", A[i][j]);
        }
        printf("\n");
        if (f) fprintf(f, "\n");
    }
}

void identity_matrix(double I[N][N], int n) {
    memset(I, 0, sizeof(double)*N*N);
    for (int i = 0; i < n; i++)
        I[i][i] = 1.0;
}

void copy_matrix(double src[N][N], double dest[N][N], int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            dest[i][j] = src[i][j];
}

// --- Gauss-Jordan ---
void inverse_gauss_jordan(double A[N][N], double inv[N][N], int n) {
    double aug[N][2*N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            aug[i][j] = A[i][j];
            aug[i][j+n] = (i == j) ? 1.0 : 0.0;
        }

    for (int i = 0; i < n; i++) {
        if (aug[i][i] == 0.0) {
            int swap = -1;
            for (int j = i+1; j < n; j++)
                if (aug[j][i] != 0.0) { swap = j; break; }
            if (swap == -1) {
                printf("Khong kha nghich (h%d = 0)\n", i);
                return;
            }
            for (int k = 0; k < 2*n; k++) {
                double temp = aug[i][k];
                aug[i][k] = aug[swap][k];
                aug[swap][k] = temp;
            }
        }

        double pivot = aug[i][i];
        for (int j = 0; j < 2*n; j++) aug[i][j] /= pivot;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            double factor = aug[j][i];
            for (int k = 0; k < 2*n; k++)
                aug[j][k] -= factor * aug[i][k];
        }
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = aug[i][j + n];
}

// --- Laplace ---
double determinant(double A[N][N], int n) {
    if (n == 1) return A[0][0];
    if (n == 2) return A[0][0]*A[1][1] - A[0][1]*A[1][0];
    double det = 0.0;
    double minor[N][N];
    for (int k = 0; k < n; k++) {
        int subi = 0;
        for (int i = 1; i < n; i++) {
            int subj = 0;
            for (int j = 0; j < n; j++) {
                if (j == k) continue;
                minor[subi][subj++] = A[i][j];
            }
            subi++;
        }
        det += (k % 2 == 0 ? 1 : -1) * A[0][k] * determinant(minor, n-1);
    }
    return det;
}

void cofactor(double A[N][N], double cof[N][N], int n) {
    double minor[N][N];
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            int subi = 0;
            for (int row = 0; row < n; row++) {
                if (row == i) continue;
                int subj = 0;
                for (int col = 0; col < n; col++) {
                    if (col == j) continue;
                    minor[subi][subj++] = A[row][col];
                }
                subi++;
            }
            cof[i][j] = pow(-1, i + j) * determinant(minor, n-1);
        }
}

void transpose(double A[N][N], double T[N][N], int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            T[i][j] = A[j][i];
}

void inverse_laplace(double A[N][N], double inv[N][N], int n) {
    double det = determinant(A, n);
    if (det == 0) {
        printf("Khong kha nghich (det = 0)\n");
        return;
    }
    double cof[N][N], adj[N][N];
    cofactor(A, cof, n);
    transpose(cof, adj, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = adj[i][j] / det;
}

// --- Newton-Schulz ---
void multiply(double A[N][N], double B[N][N], double C[N][N], int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            C[i][j] = 0;
            for (int k = 0; k < n; k++)
                C[i][j] += A[i][k] * B[k][j];
        }
}

void subtract(double A[N][N], double B[N][N], double C[N][N], int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            C[i][j] = A[i][j] - B[i][j];
}

void scalar_multiply(double A[N][N], double k, double B[N][N], int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            B[i][j] = A[i][j] * k;
}

void inverse_newton_schulz(double A[N][N], double inv[N][N], int n) {
    double normA = 0;
    for (int i = 0; i < n; i++) {
        double rowsum = 0;
        for (int j = 0; j < n; j++)
            rowsum += fabs(A[i][j]);
        if (rowsum > normA)
            normA = rowsum;
    }

    double X[N][N], I[N][N], AX[N][N], R[N][N], temp[N][N];
    identity_matrix(I, n);
    scalar_multiply(A, 1.0 / normA / normA, X, n);

    for (int iter = 0; iter < MAX_ITER; iter++) {
        multiply(A, X, AX, n);
        subtract(I, AX, R, n);
        multiply(R, X, temp, n);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                X[i][j] += temp[i][j];

        double max_diff = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (fabs(R[i][j]) > max_diff)
                    max_diff = fabs(R[i][j]);
        if (max_diff < EPSILON)
            break;
    }

    copy_matrix(X, inv, n);
}

// --- MAIN ---
int main() {
    int n, method;
    double A[N][N], result[N][N];

    printf("Nhap cap ma tran vuong n (toi da %d): ", N);
    scanf("%d", &n);
    printf("Nhap ma tran %dx%d:\n", n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            scanf("%lf", &A[i][j]);

    printf("\nChon phuong phap tinh nghich dao:\n");
    printf("1. Gauss-Jordan\n2. Laplace\n3. Newton-Schulz\n");
    printf("Lua chon cua ban: ");
    scanf("%d", &method);

    FILE *f = fopen("inverse_matrix.txt", "w");
    if (!f) {
        printf("Khong the mo tep de ghi ket qua.\n");
        return 1;
    }

    switch (method) {
        case 1:
            printf("\n--- Nghich dao (Gauss-Jordan) ---\n");
            fprintf(f, "--- Nghich dao (Gauss-Jordan) ---\n");
            inverse_gauss_jordan(A, result, n);
            break;
        case 2:
            printf("\n--- Nghich dao (Laplace) ---\n");
            fprintf(f, "--- Nghich dao (Laplace) ---\n");
            inverse_laplace(A, result, n);
            break;
        case 3:
            printf("\n--- Nghich dao (Newton-Schulz) ---\n");
            fprintf(f, "--- Nghich dao (Newton-Schulz) ---\n");
            inverse_newton_schulz(A, result, n);
            break;
        default:
            printf("Lua chon khong hop le.\n");
            fclose(f);
            return 1;
    }

    print_matrix(result, n, f);
    fclose(f);
    printf("\nDa luu ket qua vao tep inverse_matrix.txt\n");

    return 0;
}
