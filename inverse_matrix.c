#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define MAX_SIZE 20
#define EPSILON 1e-10

void intro(){
    printf("|==========================================================|\n");
    printf("|                 DO AN LAP TRINH TINH TOAN                |\n");
    printf("|                                                          |\n");
    printf("|                  TINH MA TRAN NGHICH DAO                 |\n");
    printf("|__________________________________________________________|\n");
    printf("|         SINH VIEN: NGUYEN TIEN DAT & NGUYEN ANH HOA      |\n");
    printf("|                                                          |\n");
    printf("|               HUONG DAN: TS.NGUYEN VAN HIEU              | \n");
    printf("|==========================================================|\n");
    printf("|                   CHUONG TRINH THUC HIEN                 | \n");
    printf("|==========================================================|\n\n");
}

void print_matrix(double **matrix, int n, const char *name) {
    printf("%s:\n", name);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%8.4f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_augmented_matrix(double **matrix, int n, const char *name) {
    printf("%s:\n", name);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 2*n; j++) {
            if (j == n) printf("| "); // Thêm dấu phân cách giữa ma trận vuông và ma trận bổ sung
            printf("%8.4f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

double determinant(double **matrix, int n) {
    if (n == 1) return matrix[0][0];
    if (n == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    
    double det = 0;
    double **minor = malloc((n-1) * sizeof(double *));
    for (int i = 0; i < n-1; i++)
        minor[i] = malloc((n-1) * sizeof(double));
    
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
        printf("Nhap hang cua ma tran (1~20): ");
        if (scanf("%d", n) != 1 || *n < 1 || *n > MAX_SIZE) {
            printf("Hang ma tran phai tu 1 den 20!\n");
            while (getchar() != '\n');
            continue;
        }
        break;
    }
    
    double **matrix = malloc(*n * sizeof(double *));
    for (int i = 0; i < *n; i++) matrix[i] = malloc(*n * sizeof(double));
    
    printf("Nhap ma tran cap %dx%d:\n", *n, *n);
    for (int i = 0; i < *n; i++) {
        printf("Nhap hang %d (cac phan tu cach nhau boi dau cach): ", i+1);
        for (int j = 0; j < *n; j++) {
            if (scanf("%lf", &matrix[i][j]) != 1) {
                printf("Vui long nhap so hop le!\n");
                while (getchar() != '\n');
                j--;
            }
        }
    }
    return matrix;
}

double **gauss_jordan(double **matrix, int n, int step_by_step) {
    double **augmented = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        augmented[i] = malloc(2*n * sizeof(double));
        for (int j = 0; j < n; j++) {
            augmented[i][j] = matrix[i][j];
            augmented[i][j+n] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    if (step_by_step) {
        printf("Bat dau phuong phap Gauss-Jordan:\n");
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
            printf("Buoc %d: Chuan hoa hang %d\n", i+1, i+1);
            print_augmented_matrix(augmented, n, "Ma tran");
        }
        
        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = augmented[k][i];
                for (int j = 0; j < 2*n; j++) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
                if (step_by_step) {
                    printf("Buoc %d.%d: Bien doi hang %d\n", i+1, k+1, k+1);
                    print_augmented_matrix(augmented, n, "Ma tran");
                }
            }
        }
    }
    
    double **inverse = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        inverse[i] = malloc(n * sizeof(double));
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
    
    double **adj = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) adj[i] = malloc(n * sizeof(double));
    
    if (step_by_step) printf("Bat dau phuong phap Laplace:\n");
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double **minor = malloc((n-1) * sizeof(double *));
            for (int k = 0; k < n-1; k++) minor[k] = malloc((n-1) * sizeof(double));
            
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
            if (step_by_step) printf("Phan tu phu hop (%d,%d): %.4f\n", i+1, j+1, adj[j][i]);
            
            for (int k = 0; k < n-1; k++) free(minor[k]);
            free(minor);
        }
    }
    
    double **inverse = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        inverse[i] = malloc(n * sizeof(double));
        for (int j = 0; j < n; j++) inverse[i][j] = adj[i][j] / det;
    }
    
    if (step_by_step) {
        print_matrix(adj, n, "Ma tran phu hop");
        printf("Dinh thuc: %.4f\n", det);
        print_matrix(inverse, n, "Ma tran nghich dao");
    }
    
    for (int i = 0; i < n; i++) free(adj[i]);
    free(adj);
    return inverse;
}

double **newton_schulz(double **matrix, int n, int step_by_step) {
    // Khởi tạo ma trận X0 = A^T / ||A||_2^2
    double **X = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        X[i] = malloc(n * sizeof(double));
    }
    
    // Tính norm Frobenius của A
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
    
    // Khởi tạo X0 = A^T / ||A||_2^2
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            X[i][j] = matrix[j][i] / (norm * norm); // Transpose của A chia cho bình phương norm
        }
    }
    
    if (step_by_step) {
        printf("Bat dau phuong phap Newton-Schulz:\n");
        print_matrix(X, n, "Ma tran khoi dau X0");
    }
    
    // Ma trận tạm để lưu kết quả trung gian
    double **temp = malloc(n * sizeof(double *));
    double **AX = malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        temp[i] = malloc(n * sizeof(double));
        AX[i] = malloc(n * sizeof(double));
    }
    
    // Vòng lặp Newton-Schulz: X_{k+1} = X_k (2I - A X_k)
    for (int iter = 0; iter < 100; iter++) {
        // Tính A * X_k
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                AX[i][j] = 0;
                for (int k = 0; k < n; k++) {
                    AX[i][j] += matrix[i][k] * X[k][j];
                }
            }
        }
        
        // Tính 2I - A X_k
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                temp[i][j] = (i == j ? 2.0 : 0.0) - AX[i][j];
            }
        }
        
        // Tính X_{k+1} = X_k * (2I - A X_k)
        double **newX = malloc(n * sizeof(double *));
        for (int i = 0; i < n; i++) {
            newX[i] = malloc(n * sizeof(double));
            for (int j = 0; j < n; j++) {
                newX[i][j] = 0;
                for (int k = 0; k < n; k++) {
                    newX[i][j] += X[i][k] * temp[k][j];
                }
            }
        }
        
        // Sao chép newX vào X
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                X[i][j] = newX[i][j];
            }
        }
        
        if (step_by_step) {
            printf("Vong lap %d:\n", iter+1);
            print_matrix(X, n, "Ma tran X");
        }
        
        // Kiểm tra hội tụ: ||A X - I|| < EPSILON
        double error = 0;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                double target = (i == j) ? 1.0 : 0.0;
                error += fabs(AX[i][j] - target);
            }
        }
        
        // Giải phóng newX
        for (int i = 0; i < n; i++) free(newX[i]);
        free(newX);
        
        if (error < EPSILON) break;
    }
    
    if (step_by_step) {
        print_matrix(X, n, "Ma tran nghich dao");
    }
    
    // Giải phóng bộ nhớ
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
        printf("Khong the mo tep de ghi!\n");
        return;
    }
    
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(file, "%.4f ", matrix[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
    printf("Da luu ket qua vao tep %s\n", filename);
}

int main() {
    intro();
    while (1) {
        int n;
        double **matrix = input_matrix(&n);
        double det = determinant(matrix, n);
        
        if (fabs(det) < EPSILON) {
            printf("Ma tran khong kha nghich vi dinh thuc bang 0!\n");
            for (int i = 0; i < n; i++) free(matrix[i]);
            free(matrix);
            continue;
        }
        
        while (1) {
            printf("\nMENU:\n");
            printf("1. Tinh ma tran nghich dao theo phuong phap Gauss-Jordan\n");
            printf("2. Tinh ma tran nghich dao theo phuong phap Laplace\n");
            printf("3. Tinh ma tran nghich dao theo phuong phap Newton-Schulz\n");
            printf("4. Thoat\n");
            printf("Lua chon cua ban: ");
            
            int choice;
            scanf("%d", &choice);
            if (choice == 4) {
                for (int i = 0; i < n; i++) free(matrix[i]);
                free(matrix);
                printf("Tam biet!\n");
                return 0;
            }
            
            if (choice < 1 || choice > 4) {
                printf("Lua chon khong hop le!\n");
                continue;
            }
            
            char step;
            printf("Ban co muon in ra tung buoc? (y/n): ");
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
                printf("Khong the tinh ma tran nghich dao!\n");
                continue;
            }
            
            if (!step_by_step) print_matrix(inverse, n, "Ma tran nghich dao");
            
            char save;
            printf("Ban co muon luu ket qua vao tep? (y/n): ");
            scanf(" %c", &save);
            if (save == 'y' || save == 'Y') save_to_file(inverse, n, method_name);
            
            for (int i = 0; i < n; i++) free(inverse[i]);
            free(inverse);
            
            char cont;
            printf("Ban co muon tiep tuc? (y/n): ");
            scanf(" %c", &cont);
            if (cont == 'n' || cont == 'N') {
                for (int i = 0; i < n; i++) free(matrix[i]);
                free(matrix);
                printf("Tam biet!\n");
                return 0;
            }
            break;
        }
    }
    return 0;
}
