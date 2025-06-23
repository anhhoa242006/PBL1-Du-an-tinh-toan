#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define COLOR_BRIGHT_CYAN "\x1b[96m"
#define COLOR_BRIGHT_YELLOW "\x1b[93m"
#define COLOR_BRIGHT_GREEN "\x1b[92m"
#define COLOR_BRIGHT_RED "\x1b[91m"
#define COLOR_BRIGHT_BLUE "\x1b[94m"
#define COLOR_BRIGHT_MAGENTA "\x1b[95m"
#define COLOR_RESET "\x1b[0m"

#define MAX_SIZE 20
#define EPSILON 1e-10

void intro() {
    printf("%s|==========================================================|\n", COLOR_BRIGHT_CYAN);
    printf("%s|                 DO AN LAP TRINH TINH TOAN                |\n", COLOR_BRIGHT_YELLOW);
    printf("%s|                                                          |\n", COLOR_BRIGHT_GREEN);
    printf("%s|                  TINH MA TRAN NGHICH DAO                 |\n", COLOR_BRIGHT_RED);
    printf("%s|__________________________________________________________|\n", COLOR_BRIGHT_BLUE);
    printf("%s|         SINH VIEN: NGUYEN TIEN DAT & NGUYEN ANH HOA      |\n", COLOR_BRIGHT_MAGENTA);
    printf("%s|                                                          |\n", COLOR_BRIGHT_CYAN);
    printf("%s|               HUONG DAN: TS.NGUYEN VAN HIEU              | \n", COLOR_BRIGHT_YELLOW);
    printf("%s|==========================================================|\n", COLOR_BRIGHT_GREEN);
    printf("%s|                   CHUONG TRINH THUC HIEN                 | \n", COLOR_BRIGHT_RED);
    printf("%s|==========================================================|\n%s", COLOR_BRIGHT_BLUE, COLOR_RESET);
}

void print_matrix(double **matrix, int n, const char *name) {
    printf("%s%s:\n", COLOR_BRIGHT_GREEN, name);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%8.4f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n%s", COLOR_RESET);
}

void print_augmented_matrix(double **matrix, int n, const char *name) {
    printf("%s%s:\n", COLOR_BRIGHT_GREEN, name);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 2*n; j++) {
            if (j == n) printf(" | ");
            printf("%8.4f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n%s", COLOR_RESET);
}

double determinant(double **matrix, int n) {
    if (n == 1) return matrix[0][0];
    if (n == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    
    double det = 0;
    double **minor = (double **)malloc((n-1) * sizeof(double *));
    for (int i = 0; i < n-1; i++)
        minor[i] = (double *)malloc((n-1) * sizeof(double));
    
    for (int j = 0; j < n; j++) {
        for (int row = 1; row < n; row++) {
            int col_minor = 0;
            for (int col = 0; col < n; col++) {
                if (col != j) {
                    minor[row-1][col_minor] = matrix[row][col];
                    col_minor++;
                }
            }
        }
        det += (j % 2 == 0 ? 1 : -1) * matrix[0][j] * determinant(minor, n-1);
    }
    
    for (int i = 0; i < n-1; i++) free(minor[i]);
    free(minor);
    return det;
}

double **input_matrix(int *n) {
    while (1) {
        printf("%s\nNhap hang cua ma tran (1~20): %s", COLOR_BRIGHT_BLUE, COLOR_RESET);
        if (scanf("%d", n) != 1 || *n < 1 || *n > MAX_SIZE) {
            printf("%sHang ma tran phai tu 1 den 20!\n%s", COLOR_BRIGHT_RED, COLOR_RESET);
            while (getchar() != '\n');
            continue;
        }
        break;
    }
    
    double **matrix = (double **)malloc(*n * sizeof(double *));
    for (int i = 0; i < *n; i++) matrix[i] = (double *)malloc(*n * sizeof(double));
    
    printf("%sNhap ma tran cap %dx%d:\n%s", COLOR_BRIGHT_BLUE, *n, *n, COLOR_RESET);
    for (int i = 0; i < *n; i++) {
        printf("%sNhap hang %d: %s", COLOR_BRIGHT_BLUE, i+1, COLOR_RESET);
        for (int j = 0; j < *n; j++) {
            if (scanf("%lf", &matrix[i][j]) != 1) {
                printf("%sVui long nhap so hop le!\n%s", COLOR_BRIGHT_RED, COLOR_RESET);
                while (getchar() != '\n');
                j--;
            }
        }
    }
    return matrix;
}

double **gauss_jordan(double **matrix, int n, int step_by_step) {
    double **augmented = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        augmented[i] = (double *)malloc(2*n * sizeof(double));
        for (int j = 0; j < n; j++) {
            augmented[i][j] = matrix[i][j];
            augmented[i][j+n] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    if (step_by_step) {
        printf("%s\nBat dau phuong phap Gauss-Jordan:\n%s", COLOR_BRIGHT_CYAN, COLOR_RESET);
        print_augmented_matrix(augmented, n, "Ma tran bo sung ban dau");
    }
    
    for (int i = 0; i < n; i++) {
        double pivot = augmented[i][i];
        if (fabs(pivot) < EPSILON) {
            for (int k = 0; k < n; k++) free(augmented[k]);
            free(augmented);
            return NULL;
        }
        for (int j = 0; j < 2*n; j++) augmented[i][j] /= pivot;
        
        if (step_by_step) {
            printf("%sBuoc %d: Chuan hoa hang %d\n%s", COLOR_BRIGHT_CYAN, i+1, i+1, COLOR_RESET);
            print_augmented_matrix(augmented, n, "Ma tran");
        }
        
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = augmented[k][i];
                for (int j = 0; j < 2*n; j++) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
                if (step_by_step) {
                    printf("%sBuoc %d.%d: Bien doi hang %d\n%s", COLOR_BRIGHT_CYAN, i+1, k+1, k+1, COLOR_RESET);
                    print_augmented_matrix(augmented, n, "Ma tran");
                }
            }
        }
    }
    
    double **inverse = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        inverse[i] = (double *)malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) inverse[i][j] = augmented[i][j+n];
    }
    
    if (step_by_step) {
        print_matrix(inverse, n, "Ma tran nghich dao");
    }
    
    for (int i = 0; i < n; i++) free(augmented[i]);
    free(augmented);
    return inverse;
}

double **laplace_inverse(double **matrix, int n, int step_by_step) {
    double det = determinant(matrix, n);
    if (fabs(det) < EPSILON) return NULL;
    
    double **adj = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) adj[i] = (double *)malloc(n * sizeof(double));
    
    if (step_by_step) printf("%s\nBat dau phuong phap Laplace:\n%s", COLOR_BRIGHT_CYAN, COLOR_RESET);
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double **minor = (double **)malloc((n-1) * sizeof(double *));
            for (int k = 0; k < n-1; k++) minor[k] = (double *)malloc((n-1) * sizeof(double));
            
            for (int row = 0, m_row = 0; row < n; row++) {
                if (row == i) continue;
                for (int col = 0, m_col = 0; col < n; col++) {
                    if (col == j) continue;
                    minor[m_row][m_col] = matrix[row][col];
                    m_col++;
                }
                m_row++;
            }
            
            double cofactor = (i+j) % 2 == 0 ? 1 : -1;
            adj[j][i] = cofactor * determinant(minor, n-1);
            if (step_by_step) printf("%sPhan tu phu hop (%d,%d): %.4f\n%s", COLOR_BRIGHT_CYAN, i+1, j+1, adj[j][i], COLOR_RESET);
            
            for (int k = 0; k < n-1; k++) free(minor[k]);
            free(minor);
        }
    }
    
    double **inverse = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        inverse[i] = (double *)malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) inverse[i][j] = adj[i][j] / det;
    }
    
    if (step_by_step) {
        print_matrix(adj, n, "Ma tran phu hop");
        printf("%sDinh thuc: %.4f\n%s", COLOR_BRIGHT_CYAN, det, COLOR_RESET);
        print_matrix(inverse, n, "Ma tran nghich dao");
    }
    
    for (int i = 0; i < n; i++) free(adj[i]);
    free(adj);
    return inverse;
}

double **newton_schulz(double **matrix, int n, int step_by_step) {
    double **X = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        X[i] = (double *)malloc(n * sizeof(double));
    }
    double norm = 0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            norm += matrix[i][j] * matrix[i][j];
        }
    }
    norm = sqrt(norm);
    if (norm < EPSILON) {
        for (int i = 0; i < n; i++) free(X[i]);
        free(X);
        return NULL;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            X[i][j] = matrix[j][i] / (norm * norm); 
        }
    }
    
    if (step_by_step) {
        printf("%s\nBat dau phuong phap Newton-Schulz:\n%s", COLOR_BRIGHT_CYAN, COLOR_RESET);
        print_matrix(X, n, "Ma tran khoi dau X0");
    }
    double **temp = (double **)malloc(n * sizeof(double *));
    double **AX = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        temp[i] = (double *)malloc(n * sizeof(double));
        AX[i] = (double *)malloc(n * sizeof(double));
    }
    for (int iter = 0; iter < 100; iter++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                AX[i][j] = 0;
                for (int k = 0; k < n; k++) {
                    AX[i][j] += matrix[i][k] * X[k][j];
                }
            }
        }
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                temp[i][j] = (i == j ? 2.0 : 0.0) - AX[i][j];
            }
        }
        
        double **newX = (double **)malloc(n * sizeof(double *));
        for (int i = 0; i < n; i++) {
            newX[i] = (double *)malloc(n * sizeof(double));
            for (int j = 0; j < n; j++) {
                newX[i][j] = 0;
                for (int k = 0; k < n; k++) {
                    newX[i][j] += X[i][k] * temp[k][j];
                }
            }
        }
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                X[i][j] = newX[i][j];
            }
        }
        
        if (step_by_step) {
            printf("%sVong lap %d:\n%s", COLOR_BRIGHT_CYAN, iter+1, COLOR_RESET);
            print_matrix(X, n, "Ma tran X");
        }
        
        double error = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double target = (i == j) ? 1.0 : 0.0;
                error += fabs(AX[i][j] - target);
            }
        }
        
        for (int i = 0; i < n; i++) free(newX[i]);
        free(newX);
        
        if (error < EPSILON) break;
    }
    
    if (step_by_step) {
        print_matrix(X, n, "Ma tran nghich dao");
    }
    
    for (int i = 0; i < n; i++) {
        free(temp[i]);
        free(AX[i]);
    }
    free(temp);
    free(AX);
    
    return X;
}

void save_to_file(double **matrix, int n, const char *method) {
    char filename[50];
    snprintf(filename, sizeof(filename), "inverse_%s.txt", method);
    FILE *file = fopen(filename, "w");
    if (!file) {
        printf("%sKhong the mo tep de ghi!\n%s", COLOR_BRIGHT_RED, COLOR_RESET);
        return;
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(file, "%.4f ", matrix[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
    printf("%sDa luu ket qua vao tep.%s\n%s", COLOR_BRIGHT_MAGENTA, filename, COLOR_RESET);
}

int main() {
    intro();
    while (1) {
        int n;
        double **matrix = input_matrix(&n);
        double det = determinant(matrix, n);
        
        if (fabs(det) < EPSILON) {
            printf("%sMa tran khong kha nghich vi dinh thuc bang 0!\n%s", COLOR_BRIGHT_RED, COLOR_RESET);
            for (int i = 0; i < n; i++) free(matrix[i]);
            free(matrix);
            continue;
        }
        
        while (1) {
            printf("%s\nMENU:\n", COLOR_BRIGHT_YELLOW);
            printf("1. Tinh ma tran nghich dao theo phuong phap Gauss-Jordan\n");
            printf("2. Tinh ma tran nghich dao theo phuong phap Laplace\n");
            printf("3. Tinh ma tran nghich dao theo phuong phap Newton-Schulz\n");
            printf("4. Thoat\n");
            printf("\nLua chon cua ban: %s", COLOR_RESET);
            
            int choice;
            scanf("%d", &choice);
            if (choice == 4) {
                for (int i = 0; i < n; i++) free(matrix[i]);
                free(matrix);
                printf("%sCam on ban da su dung chuong trinh.\n%s", COLOR_BRIGHT_CYAN, COLOR_RESET);
                return 0;
            }
            
            if (choice < 1 || choice > 4) {
                printf("%sLua chon khong hop le!\n%s", COLOR_BRIGHT_RED, COLOR_RESET);
                continue;
            }
            
            char step;
            printf("%sBan co muon in ra tung buoc? (y/n): %s", COLOR_BRIGHT_BLUE, COLOR_RESET);
            scanf(" %c", &step);
            int step_by_step = (step == 'y' || step == 'Y');
            
            double **inverse = NULL;
            const char *method_name = "";
            
            switch (choice) {
                case 1:
                    inverse = gauss_jordan(matrix, n, step_by_step);
                    method_name = "Gauss-Jordan";
                    break;
                case 2:
                    inverse = laplace_inverse(matrix, n, step_by_step);
                    method_name = "Laplace";
                    break;
                case 3:
                    inverse = newton_schulz(matrix, n, step_by_step);
                    method_name = "Newton-Schulz";
                    break;
            }
            
            if (!inverse) {
                printf("%sKhong the tinh ma tran nghich dao!\n%s", COLOR_BRIGHT_RED, COLOR_RESET);
                continue;
            }
            
            if (!step_by_step) print_matrix(inverse, n, "Ma tran nghich dao");
            
            char save;
            printf("%sBan co muon luu ket qua vao tep? (y/n): %s", COLOR_BRIGHT_BLUE, COLOR_RESET);
            scanf(" %c", &save);
            if (save == 'y' || save == 'Y') save_to_file(inverse, n, method_name);
            
            for (int i = 0; i < n; i++) free(inverse[i]);
            free(inverse);
            
            char cont;
            printf("%s\nBan co muon tiep tuc? (y/n): %s", COLOR_BRIGHT_BLUE, COLOR_RESET);
            scanf(" %c", &cont);
            if (cont == 'n' || cont == 'N') {
                for (int i = 0; i < n; i++) free(matrix[i]);
                free(matrix);
                printf("%s\nCam on ban da su dung chuong trinh.\n%s", COLOR_BRIGHT_CYAN, COLOR_RESET);
                return 0;
            }
            break;
        }
    }
    return 0;
}
