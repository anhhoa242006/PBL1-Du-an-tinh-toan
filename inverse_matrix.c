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
    if (!matrix) {
        printf("Loi cap phat bo nho cho ma tran.\n");
        exit(1);
    }
    for (int i = 0; i < n; i++) {
        matrix[i] = (float *)calloc(n, sizeof(float));
        if (!matrix[i]) {
            printf("Loi cap phat bo nho cho hang %d.\n", i);
            exit(1);
        }
    }
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
        for (int j = 0; j < n; j++) {
            while (scanf("%f", &a[i][j]) != 1) {
                printf("Gia tri khong hop le, vui long nhap lai: ");
                while (getchar() != '\n'); // Clear input buffer
            }
        }
}

void printMatrix(float **matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fabs(matrix[i][j]) < EPSILON)
                printf("%8.3f ", 0.0);
            else
                printf("%8.3f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void printAugmentedMatrix(float **a, float **inv, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
            printf("%8.3f ", fabs(a[i][j]) < EPSILON ? 0.0 : a[i][j]);
        printf(" | ");
        for (int j = 0; j < n; j++)
            printf("%8.3f ", fabs(inv[i][j]) < EPSILON ? 0.0 : inv[i][j]);
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
        // Partial pivoting: find row with largest pivot
        int maxRow = i;
        float maxVal = fabs(tmp[i][i]);
        for (int k = i + 1; k < n; k++) {
            if (fabs(tmp[k][i]) > maxVal) {
                maxVal = fabs(tmp[k][i]);
                maxRow = k;
            }
        }
        if (fabs(maxVal) < EPSILON) {
            printf("Khong the nghich dao (cot %d toan 0).\n", i);
            freeMatrix(tmp, n);
            return 0;
        }
        if (maxRow != i) {
            float *t = tmp[i]; tmp[i] = tmp[maxRow]; tmp[maxRow] = t;
            t = inv[i]; inv[i] = inv[maxRow]; inv[maxRow] = t;
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
            printAugmentedMatrix(tmp, inv, n);
            printf("\n");
        }
    }

    freeMatrix(tmp, n);
    return 1;
}

float determinant(float **a, int n, int printSteps) {
    if (n == 1) return a[0][0];
    if (n == 2) return a[0][0] * a[1][1] - a[0][1] * a[1][0];

    float det = 0;
    float **sub = createMatrix(n - 1);

    for (int k = 0; k < n; k++) {
        if (fabs(a[0][k]) < EPSILON) continue; // Skip near-zero elements
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
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            int row = 0;
            for (int ii = 0; ii < n; ii++) {
                if (ii == i) continue;
                int col = 0;
                for (int jj = 0; jj < n; jj++) {
                    if (jj == j) continue;
                    sub[row][col++] = a[ii][jj];
                }
                row++;
            }

            float cof = (((i + j) % 2 == 0) ? 1 : -1) * determinant(sub, n - 1, 0);
            adj[i][j] = cof;

            if (printSteps) {
                printf("Minor (%d,%d):\n", i, j);
                printMatrix(sub, n - 1);
                printf("Cofactor: %.3f\n\n", cof);
            }
        }
    }
    freeMatrix(sub, n - 1);
}

int inverseLaplace(float **a, float **inv, int n, int printSteps) {
    float det = determinant(a, n, printSteps);
    if (fabs(det) < EPSILON) {
        printf("Ma tran khong kha nghich (dinh thuc gan bang 0).\n");
        return 0;
    }

    float **adj = createMatrix(n);
    adjugate(a, adj, n, printSteps);

    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            inv[i][j] = adj[j][i] / det;

    freeMatrix(adj, n);
    return 1;
}

int inverseNewtonSchulz(float **a, float **inv, int n, int printSteps) {
    float **X = createMatrix(n);
    float **I = createMatrix(n);
    float **AX = createMatrix(n);
    float **R = createMatrix(n);
    float **newX = createMatrix(n);

    for (int i = 0; i < n; i++) I[i][i] = 1.0f;

    // Compute Frobenius norm for initial guess
    float norm = 0.0f;
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            norm += a[i][j] * a[i][j];
    norm = sqrt(norm);
    float scale = (norm != 0) ? 1.0f / (norm * norm) : 1.0f;

    // Initial guess: X0 = A^T / ||A||_F^2
    for (int i = 0; i < n; i++)
        for (int j = 0; j < n; j++)
            X[i][j] = a[j][i] * scale;

    int converged = 0;
    for (int iter = 0; iter < MAX_ITER; iter++) {
        // Compute AX
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                AX[i][j] = 0;
                for (int k = 0; k < n; k++)
                    AX[i][j] += a[i][k] * X[k][j];
            }

        // R = 2I - AX
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                R[i][j] = 2.0f * I[i][j] - AX[i][j];

        // newX = X * R
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                newX[i][j] = 0;
                for (int k = 0; k < n; k++)
                    newX[i][j] += X[i][k] * R[k][j];
            }

        // Check convergence: ||I - A*X||_F
        float error = 0;
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++) {
                float val = (i == j) ? 1.0f : 0.0f;
                val -= AX[i][j];
                error += val * val;
            }
        error = sqrt(error);

        if (printSteps) {
            printf("Lan lap %d (error: %.6f):\n", iter + 1, error);
            printMatrix(newX, n);
        }
        
        if (error < EPSILON) {
            converged = 1;
            break;
        }

        copyMatrix(newX, X, n);
    }

    if (!converged) {
        printf("Phuong phap Newton-Schulz khong hoi tu sau %d lan lap.\n", MAX_ITER);
        freeMatrix(X, n); freeMatrix(I, n);
        freeMatrix(AX, n); freeMatrix(R, n);
        freeMatrix(newX, n);
        return 0;
    }

    copyMatrix(X, inv, n);

    freeMatrix(X, n); freeMatrix(I, n);
    freeMatrix(AX, n); freeMatrix(R, n);
    freeMatrix(newX, n);
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
            fprintf(f, "%8.3f ", fabs(a[i][j]) < EPSILON ? 0.0 : a[i][j]);
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
            while (scanf("%d", &n) != 1 || n <= 0) {
                printf("Cap ma tran khong hop le, vui long nhap lai: ");
                while (getchar() != '\n');
            }
        } while (n <= 0);

        float **A = createMatrix(n);
        float **inv = createMatrix(n);
        inputMatrix(A, n);

        do {
            printf("\nChon phuong phap tinh nghich dao:\n1. Gauss-Jordan\n2. Laplace\n3. Newton-Schulz\nLua chon cua ban (1-3): ");
            while (scanf("%d", &method) != 1 || method < 1 || method > 3) {
                printf("Lua chon khong hop le, vui long nhap lai: ");
                while (getchar() != '\n');
            }
        } while (method < 1 || method > 3);

        printf("Hien thi tung buoc tinh toan? (0: Khong, 1: Co): ");
        while (scanf("%d", &printSteps) != 1 || (printSteps != 0 && printSteps != 1)) {
            printf("Lua chon khong hop le, vui long nhap lai: ");
            while (getchar() != '\n');
        }

        int success = 0;
        switch (method) {
            case 1: success = inverseGaussJordan(A, inv, n, printSteps); strcpy(methodName, "Gauss-Jordan"); break;
            case 2: success = inverseLaplace(A, inv, n, printSteps); strcpy(methodName, "Laplace"); break;
            case 3: success = inverseNewtonSchulz(A, inv, n, printSteps); strcpy(methodName, "Newton-Schulz"); break;
        }

        if (success) {
            printf("Ma tran nghich dao la:\n");
            printMatrix(inv, n);
            printf("Ban co muon luu ket qua vao tep? (0: Khong, 1: Co): ");
            while (scanf("%d", &save) != 1 || (save != 0 && save != 1)) {
                printf("Lua chon khong hop le, vui long nhap lai: ");
                while (getchar() != '\n');
            }
            if (save) saveToFile(inv, n, "inverse_matrix.txt", methodName);
        }

        freeMatrix(A, n);
        freeMatrix(inv, n);

        printf("Ban co muon tiep tuc chuong trinh? (0: Dung, 1: Tiep tuc): ");
        while (scanf("%d", &again) != 1 || (again != 0 && again != 1)) {
            printf("Lua chon khong hop le, vui long nhap lai: ");
            while (getchar() != '\n');
        }
    } while (again);

    printf("\nCam on ban da su dung chuong trinh.\n");
    return 0;
}
