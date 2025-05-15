#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#define EPSILON 1e-6
#define MAX_ITER 100

void intro() {
    printf("***********************************************************\n");
    printf("                 DO AN LAP TRINH TINH TOAN                 \n");
    printf("***********************************************************\n");
    printf("         SINH VIEN: NGUYEN TIEN DAT & NGUYEN ANH HOA       \n");
    printf("         HUONG DAN: TS.NGUYEN VAN HIEU                     \n");
    printf("===========================================================\n");
    printf("                   CHUONG TRINH THUC HIEN                  \n");
    printf("===========================================================\n\n");
}

float **createMatrix(int n) {
    float **matrix = (float **)calloc(n, sizeof(float *));
    for (int i = 0; i < n; i++)
        matrix[i] = (float *)calloc(n, sizeof(float));
    return matrix;
}

void freeMatrix(float **matrix, int n) {
    for (int i = 0; i < n; i++) free(matrix[i]);
    free(matrix);
}

void copyMatrix(float **src, float **dest, int n) {
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            dest[i][j] = src[i][j];
}

void inputMatrix(float **a, int n) {
    printf("Nhap ma tran %dx%d:\n", n, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            scanf("%f", &a[i][j]);
}

void printMatrix(float **a, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%8.3f ", a[i][j]);
        printf("\n");
    }
}

int inverseGaussJordan(float **a, float **inv, int n, int printSteps) {
    float **tmp = createMatrix(n);
    copyMatrix(a, tmp, n);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = (i == j) ? 1.0f : 0.0f;

    for (int i = 0; i < n; i++) {
        if (fabs(tmp[i][i]) < EPSILON) {
            printf("Khong the nghich dao (pivot = 0)\n");
            freeMatrix(tmp, n);
            return 0;
        }

        float pivot = tmp[i][i];
        for (int j = 0; j < n; j++) {
            tmp[i][j] /= pivot;
            inv[i][j] /= pivot;
        }

        for (int k = 0; k < n; k++) {
            if (k == i) continue;
            float factor = tmp[k][i];
            for (int j = 0; j < n; j++) {
                tmp[k][j] -= factor * tmp[i][j];
                inv[k][j] -= factor * inv[i][j];
            }
        }

        if (printSteps) {
            printf("Buoc %d:\n", i + 1);
            printMatrix(inv, n);
        }
    }

    freeMatrix(tmp, n);
    return 1;
}

float determinant(float **a, int n, int printSteps) {
    if (n == 1) return a[0][0];
    if (n == 2) return a[0][0]*a[1][1] - a[0][1]*a[1][0];

    float det = 0;
    float **sub = createMatrix(n - 1);

    for (int k = 0; k < n; k++) {
        for (int i = 1; i < n; i++) {
            int col = 0;
            for (int j = 0; j < n; j++)
                if (j != k)
                    sub[i-1][col++] = a[i][j];
        }
        float cofactor = ((k % 2 == 0) ? 1 : -1) * a[0][k] * determinant(sub, n - 1, printSteps);
        det += cofactor;

        if (printSteps) {
            printf("Minor tai cot %d:\n", k);
            printMatrix(sub, n - 1);
            printf("Cofactor: %.3f\n", cofactor);
        }
    }

    freeMatrix(sub, n - 1);
    return det;
}

void adjugate(float **a, float **adj, int n, int printSteps) {
    float **sub = createMatrix(n - 1);
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++) {
            int row = 0, col;
            for (int ii = 0; ii < n; ii++) {
                if (ii == i) continue;
                col = 0;
                for (int jj = 0; jj < n; jj++) {
                    if (jj == j) continue;
                    sub[row][col++] = a[ii][jj];
                }
                row++;
            }
            float cof = (((i + j) % 2 == 0) ? 1 : -1) * determinant(sub, n - 1, printSteps);
            adj[j][i] = cof;

            if (printSteps) {
                printf("Minor (%d,%d):\n", i, j);
                printMatrix(sub, n - 1);
                printf("Cofactor: %.3f\n\n", cof);
            }
        }
    freeMatrix(sub, n - 1);
}

int inverseLaplace(float **a, float **inv, int n, int printSteps) {
    float det = determinant(a, n, printSteps);
    if (fabs(det) < EPSILON) {
        printf("Khong the nghich dao (dinh thuc = 0)\n");
        return 0;
    }
    float **adj = createMatrix(n);
    adjugate(a, adj, n, printSteps);

    if (printSteps) {
        printf("Ma tran phu hop (adjugate):\n");
        printMatrix(adj, n);
        printf("Dinh thuc: %.6f\n", det);
    }

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = adj[i][j] / det;
    freeMatrix(adj, n);
    return 1;
}

int inverseNewtonSchulz(float **a, float **inv, int n, int printSteps) {
    float **X = createMatrix(n);
    float **I = createMatrix(n);
    float **AX = createMatrix(n);
    float **T = createMatrix(n);

    float norm = 0;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            norm += a[i][j] * a[i][j];

    float scale = 1.0f / norm;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            X[i][j] = a[j][i] * scale;

    for (int i = 0; i < n; i++) I[i][i] = 1.0f;

    for (int iter = 0; iter < MAX_ITER; iter++) {
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                AX[i][j] = 0;
                for (int k = 0; k < n; k++) AX[i][j] += a[i][k] * X[k][j];
                T[i][j] = I[i][j] - AX[i][j];
            }

        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                float sum = 0;
                for (int k = 0; k < n; k++) sum += X[i][k] * T[k][j];
                inv[i][j] = X[i][j] + sum;
            }

        if (printSteps) {
            printf("Lan lap %d:\n", iter + 1);
            printMatrix(inv, n);
        }

        copyMatrix(inv, X, n);
    }

    freeMatrix(X, n); freeMatrix(I, n);
    freeMatrix(AX, n); freeMatrix(T, n);
    return 1;
}

void saveToFile(float **a, int n, const char *filename, const char *methodName) {
    FILE *f = fopen(filename, "a");
    if (!f) {
        printf("Khong the mo file de ghi.\n");
        return;
    }

    time_t now = time(NULL);
    struct tm *t = localtime(&now);
    fprintf(f, "\nThoi gian: %02d-%02d-%04d %02d:%02d:%02d\n",
            t->tm_mday, t->tm_mon + 1, t->tm_year + 1900,
            t->tm_hour, t->tm_min, t->tm_sec);
    fprintf(f, "Phuong phap: %s\n", methodName);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            fprintf(f, "%8.3f ", a[i][j]);
        fprintf(f, "\n");
    }
    fprintf(f, "\n");

    fclose(f);
    printf("Da luu vao file '%s'\n", filename);
}

int main() {
    int n, method, printSteps, save, again;
    char methodName[50];
    intro();

    FILE *fout = fopen("inverse_matrix.txt", "w");
    if (!fout) {
        printf("Khong the mo file.\n");
        return 1;
    }
    fclose(fout);

    do {
        do {
            printf("Nhap cap ma tran (n > 0): ");
            scanf("%d", &n);
        } while (n <= 0);

        float **A = createMatrix(n);
        float **inv = createMatrix(n);
        inputMatrix(A, n);

        do {
            printf("\nChon phuong phap:\n1. Gauss-Jordan\n2. Laplace\n3. Newton-Schulz\nLua chon (1-3): ");
            scanf("%d", &method);
        } while (method < 1 || method > 3);

        printf("In ket qua tung buoc? (0: Khong, 1: Co): ");
        scanf("%d", &printSteps);

        int success = 0;
        switch (method) {
            case 1: success = inverseGaussJordan(A, inv, n, printSteps); strcpy(methodName, "Gauss-Jordan"); break;
            case 2: success = inverseLaplace(A, inv, n, printSteps); strcpy(methodName, "Laplace"); break;
            case 3: success = inverseNewtonSchulz(A, inv, n, printSteps); strcpy(methodName, "Newton-Schulz"); break;
        }

        if (success) {
            printf("Ma tran nghich dao:\n");
            printMatrix(inv, n);
            printf("Luu vao tep? (0: Khong, 1: Co): ");
            scanf("%d", &save);
            if (save) saveToFile(inv, n, "inverse_matrix.txt", methodName);
        }

        freeMatrix(A, n);
        freeMatrix(inv, n);

        printf("Tiep tuc? (0: Dung, 1: Tiep tuc): ");
        scanf("%d", &again);
    } while (again);

    return 0;
}
