#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <windows.h>

#define COLOR_BRIGHT_CYAN "\x1b[96m"   
#define COLOR_BRIGHT_YELLOW "\x1b[93m" 
#define COLOR_BRIGHT_RED "\x1b[91m"    
#define COLOR_BRIGHT_GREEN "\x1b[92m"  
#define COLOR_BRIGHT_MAGENTA "\x1b[95m" 
#define COLOR_RESET "\x1b[0m"

#define MAX_SIZE 20
#define EPSILON 1e-10

// Viền đơn giản
void print_border(const char *color, int length) {
    for (int i = 0; i < length; i++) printf("%s-%s", color, (i == length - 1) ? "\n" : "");
    printf("%s", COLOR_RESET);
}

// Tiêu đề căn giữa
void print_title(const char *color, const char *title, int length) {
    int len = strlen(title);
    int padding = (length - len) / 2;
    printf("%s%*s%s%*s\n%s", color, padding, "", title, padding + (len % 2), "", COLOR_RESET);
}

// Hàm gioi thieu
void intro() {
    printf("%s|==========================================================|\n", COLOR_BRIGHT_CYAN);
    printf("%s|                 DO AN LAP TRINH TINH TOAN                |\n", COLOR_BRIGHT_YELLOW);
    printf("%s|                                                          |\n", COLOR_BRIGHT_CYAN);
    printf("%s|                  TINH MA TRAN NGHICH DAO                 |\n", COLOR_BRIGHT_RED);
    printf("%s|__________________________________________________________|\n", COLOR_BRIGHT_CYAN);
    printf("%s|         SINH VIEN: NGUYEN TIEN DAT & NGUYEN ANH HOA      |\n", COLOR_BRIGHT_YELLOW);
    printf("%s|                                                          |\n", COLOR_BRIGHT_CYAN);
    printf("%s|               HUONG DAN: TS.NGUYEN VAN HIEU              | \n", COLOR_BRIGHT_YELLOW);
    printf("%s|==========================================================|\n", COLOR_BRIGHT_CYAN);
    printf("%s|                   CHUONG TRINH THUC HIEN                 | \n", COLOR_BRIGHT_YELLOW);
    printf("%s|==========================================================|\n%s", COLOR_BRIGHT_CYAN, COLOR_RESET);
}

// In ma trận với căn lề thu hẹp
void print_matrix(double **matrix, int n, const char *name) {
    int max_len = 50;
    print_border(COLOR_BRIGHT_CYAN, max_len);
    print_title(COLOR_BRIGHT_CYAN, name, max_len);
    for (int i = 0; i < n; i++) {
        printf("%s     ", COLOR_BRIGHT_YELLOW);
        for (int j = 0; j < n; j++) {
            printf("%7.2f", matrix[i][j]); 
        }
        printf("\n");
    }
    print_border(COLOR_BRIGHT_CYAN, max_len);
}

// In ma trận bổ sung với căn lề thu hẹp
void print_augmented_matrix(double **matrix, int n, const char *name) {
    int max_len = 50;
    print_border(COLOR_BRIGHT_CYAN, max_len);
    print_title(COLOR_BRIGHT_CYAN, name, max_len);
    for (int i = 0; i < n; i++) {
        printf("%s     ", COLOR_BRIGHT_YELLOW);
        for (int j = 0; j < 2*n; j++) {
            if (j == n) printf(" |");
            printf("%7.2f", matrix[i][j]); 
        }
        printf("\n");
    }
    print_border(COLOR_BRIGHT_CYAN, max_len);
}

// Tính định thức
double determinant(double **matrix, int n) {
    if (n == 1) return matrix[0][0];
    if (n == 2) return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];

    double det = 0;
    double **minor = (double **)malloc((n-1) * sizeof(double *));
    for (int i = 0; i < n-1; i++)
        minor[i] = (double *)malloc((n-1) * sizeof(double));

    for (int j = 0; j < n; j++) {
        for (int row = 0, m_row = 0; row < n; row++) {
            if (row == 0) continue;
            for (int col = 0, m_col = 0; col < n; col++) {
                if (col != j) {
                    minor[m_row][m_col] = matrix[row][col];
                    m_col++;
                }
            }
            m_row++;
        }
        det += (j % 2 == 0 ? 1 : -1) * matrix[0][j] * determinant(minor, n-1);
    }

    for (int i = 0; i < n-1; i++) free(minor[i]);
    free(minor);
    return det;
}

// Nhập ma trận
double **input_matrix(int *n) {
    int max_len = 50;
    print_border(COLOR_BRIGHT_CYAN, max_len);
    print_title(COLOR_BRIGHT_CYAN, "Nhap Ma Tran", max_len);
    printf("%sCap (1-%d): ", COLOR_BRIGHT_YELLOW, MAX_SIZE);
    while (1) {
        if (scanf("%d", n) != 1 || *n < 1 || *n > MAX_SIZE) {
            printf("%sLoi: Cap tu 1 den %d!\n", COLOR_BRIGHT_RED, MAX_SIZE);
            printf("%sNhap lai: ", COLOR_BRIGHT_YELLOW);
            while (getchar() != '\n');
            continue;
        }
        while (getchar() != '\n');
        break;
    }

    double **matrix = (double **)malloc(*n * sizeof(double *));
    for (int i = 0; i < *n; i++) matrix[i] = (double *)malloc(*n * sizeof(double));

    printf("%sMa tran %dx%d (vi du: 1.5 -2.0):\n", COLOR_BRIGHT_YELLOW, *n, *n);
    for (int i = 0; i < *n; i++) {
        printf("%sHang %d: ", COLOR_BRIGHT_YELLOW, i+1);
        for (int j = 0; j < *n; j++) {
            if (scanf("%lf", &matrix[i][j]) != 1) {
                printf("%sLoi: Nhap so thuc!\n", COLOR_BRIGHT_RED);
                printf("%sHang %d: ", COLOR_BRIGHT_YELLOW, i+1);
                while (getchar() != '\n');
                j--;
            }
        }
        while (getchar() != '\n');
    }
    print_border(COLOR_BRIGHT_CYAN, max_len);
    return matrix;
}

// Gauss-Jordan
double **gauss_jordan(double **matrix, int n, int step_by_step, double *exec_time) {
    clock_t start = clock();
    double **augmented = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) {
        augmented[i] = (double *)malloc(2*n * sizeof(double));
        for (int j = 0; j < n; j++) {
            augmented[i][j] = matrix[i][j];
            augmented[i][j+n] = (i == j) ? 1.0 : 0.0;
        }
    }

    if (step_by_step) {
        print_augmented_matrix(augmented, n, "Ma tran ban dau");
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
            char step_msg[50];
            snprintf(step_msg, sizeof(step_msg), "Buoc %d: Hang %d", i+1, i+1);
            print_border(COLOR_BRIGHT_CYAN, 50);
            print_title(COLOR_BRIGHT_CYAN, step_msg, 50);
            print_augmented_matrix(augmented, n, "Ma tran");
        }

        for (int k = 0; k < n; k++) {
            if (k != i) {
                double factor = augmented[k][i];
                for (int j = 0; j < 2*n; j++) {
                    augmented[k][j] -= factor * augmented[i][j];
                }
                if (step_by_step) {
                    char step_msg[50];
                    snprintf(step_msg, sizeof(step_msg), "Buoc %d.%d: Hang %d", i+1, k+1, k+1);
                    print_border(COLOR_BRIGHT_CYAN, 50);
                    print_title(COLOR_BRIGHT_CYAN, step_msg, 50);
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

    *exec_time = ((double)(clock() - start)) / CLOCKS_PER_SEC;
    return inverse;
}

// Laplace
double **laplace_inverse(double **matrix, int n, int step_by_step, double *exec_time) {
    clock_t start = clock();
    double det = determinant(matrix, n);
    if (fabs(det) < EPSILON) return NULL;

    double **adj = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++) adj[i] = (double *)malloc(n * sizeof(double));

    if (step_by_step) {
        print_border(COLOR_BRIGHT_CYAN, 50);
        print_title(COLOR_BRIGHT_CYAN, "Laplace", 50);
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double **minor = (double **)malloc((n-1) * sizeof(double *));
            for (int k = 0; k < n-1; k++) minor[k] = (double *)malloc((n-1) * sizeof(double));

            for (int row = 0, m_row = 0; row < n; row++) {
                if (row == i) continue;
                for (int col = 0, m_col = 0; col < n; col++) {
                    if (col != j) {
                        minor[m_row][m_col] = matrix[row][col];
                        m_col++;
                    }
                }
                m_row++;
            }

            double cofactor = (i+j) % 2 == 0 ? 1 : -1;
            adj[j][i] = cofactor * determinant(minor, n-1);
            if (step_by_step) {
                printf("%sPhu hop (%d,%d): %.2f\n", COLOR_BRIGHT_CYAN, i+1, j+1, adj[j][i]);
            }

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
        printf("%sDinh thuc: %.2f\n", COLOR_BRIGHT_CYAN, det);
        print_matrix(inverse, n, "Ma tran nghich dao");
    }

    for (int i = 0; i < n; i++) free(adj[i]);
    free(adj);

    *exec_time = ((double)(clock() - start)) / CLOCKS_PER_SEC;
    return inverse;
}

// Newton-Schulz
double **newton_schulz(double **matrix, int n, int step_by_step, double *exec_time) {
    clock_t start = clock();
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
        print_border(COLOR_BRIGHT_CYAN, 50);
        print_title(COLOR_BRIGHT_CYAN, "Newton-Schulz", 50);
        print_matrix(X, n, "Ma tran X0");
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
            char iter_msg[30];
            snprintf(iter_msg, sizeof(iter_msg), "Vong lap %d", iter+1);
            print_border(COLOR_BRIGHT_CYAN, 50);
            print_title(COLOR_BRIGHT_CYAN, iter_msg, 50);
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

    *exec_time = ((double)(clock() - start)) / CLOCKS_PER_SEC;
    return X;
}

// Lưu kết quả
void save_to_file(double **matrix, int n, const char *method, double exec_time) {
    char filename[50];
    snprintf(filename, sizeof(filename), "inverse_%s.txt", method);
    FILE *file = fopen(filename, "w");
    if (!file) {
        print_border(COLOR_BRIGHT_RED, 50);
        print_title(COLOR_BRIGHT_RED, "Loi tep!", 50);
        return;
    }

    fprintf(file, "Ma tran nghich dao (%s)\n", method);
    fprintf(file, "Thoi gian: %.6f giay\n", exec_time);
    fprintf(file, "--------------------------------------------------\n");
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            fprintf(file, "%10.4f ", matrix[i][j]);
        }
        fprintf(file, "\n");
    }
    fclose(file);
    char success_msg[100];
    snprintf(success_msg, sizeof(success_msg), "Luu thanh cong vao %s", filename);
    print_border(COLOR_BRIGHT_GREEN, 50);
    print_title(COLOR_BRIGHT_GREEN, success_msg, 50);
}

// In menu với màu magenta
void print_menu() {
    int max_len = 50;
    print_border(COLOR_BRIGHT_MAGENTA, max_len);
    print_title(COLOR_BRIGHT_MAGENTA, "MENU", max_len);
    printf("%s 1. Gauss-Jordan\n", COLOR_BRIGHT_YELLOW);
    printf("%s 2. Laplace\n", COLOR_BRIGHT_YELLOW);
    printf("%s 3. Newton-Schulz\n", COLOR_BRIGHT_YELLOW);
    printf("%s 4. Thoat\n", COLOR_BRIGHT_YELLOW);
    print_border(COLOR_BRIGHT_MAGENTA, max_len);
    printf("%sLua chon cua ban (1-4): %s", COLOR_BRIGHT_YELLOW, COLOR_RESET);
}

int main() {
    // Kích hoạt ANSI
    HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
    DWORD consoleMode;
    GetConsoleMode(hConsole, &consoleMode);
    SetConsoleMode(hConsole, consoleMode | ENABLE_VIRTUAL_TERMINAL_PROCESSING);

    system("cls");
    intro();
    while (1) {
        int n;
        double **matrix = input_matrix(&n);
        double det = determinant(matrix, n);

        if (fabs(det) < EPSILON) {
            print_border(COLOR_BRIGHT_RED, 50);
            print_title(COLOR_BRIGHT_RED, "Loi: Khong kha nghich!", 50);
            for (int i = 0; i < n; i++) free(matrix[i]);
            free(matrix);
            continue;
        }

        while (1) {
            system("cls");
            print_menu();
            int choice;
            if (scanf("%d", &choice) != 1) {
                printf("%sLoi: Nhap so tu 1-4!\n%s", COLOR_BRIGHT_RED, COLOR_RESET);
                while (getchar() != '\n');
                continue;
            }
            while (getchar() != '\n');
            if (choice == 4) {
                for (int i = 0; i < n; i++) free(matrix[i]);
                free(matrix);
                system("cls");
                print_border(COLOR_BRIGHT_CYAN, 50);
                print_title(COLOR_BRIGHT_CYAN, "Cam on da su dung chuong trinh!", 50);
                print_border(COLOR_BRIGHT_CYAN, 50);
                return 0;
            }

            if (choice < 1 || choice > 4) {
                print_border(COLOR_BRIGHT_RED, 50);
                print_title(COLOR_BRIGHT_RED, "Loi: Chon sai!", 50);
                continue;
            }

            char step;
            print_border(COLOR_BRIGHT_CYAN, 50);
            printf("%sBan co muon in ra tung buoc? (y/n): ", COLOR_BRIGHT_YELLOW);
            scanf(" %c", &step);
            while (getchar() != '\n');
            int step_by_step = (step == 'y' || step == 'Y');
            print_border(COLOR_BRIGHT_CYAN, 50);

            double **inverse = NULL;
            const char *method_name = "";
            double exec_time = 0.0;

            switch (choice) {
                case 1:
                    inverse = gauss_jordan(matrix, n, step_by_step, &exec_time);
                    method_name = "Gauss-Jordan";
                    break;
                case 2:
                    inverse = laplace_inverse(matrix, n, step_by_step, &exec_time);
                    method_name = "Laplace";
                    break;
                case 3:
                    inverse = newton_schulz(matrix, n, step_by_step, &exec_time);
                    method_name = "Newton-Schulz";
                    break;
            }

            if (!inverse) {
                print_border(COLOR_BRIGHT_RED, 50);
                print_title(COLOR_BRIGHT_RED, "Khong the tinh ma tran nghich dao!", 50);
                continue;
            }

            if (!step_by_step) {
                print_matrix(inverse, n, "Ma tran nghich dao");
                char time_msg[30];
                snprintf(time_msg, sizeof(time_msg), "Thoi gian: %.4fs", exec_time);
                print_border(COLOR_BRIGHT_CYAN, 50);
                print_title(COLOR_BRIGHT_CYAN, time_msg, 50);
                print_border(COLOR_BRIGHT_CYAN, 50);
            }

            char save;
            print_border(COLOR_BRIGHT_CYAN, 50);
            printf("%sBan co muon luu ket qua vao tep? (y/n): ", COLOR_BRIGHT_YELLOW);
            scanf(" %c", &save);
            while (getchar() != '\n');
            print_border(COLOR_BRIGHT_CYAN, 50);
            if (save == 'y' || save == 'Y') save_to_file(inverse, n, method_name, exec_time);

            for (int i = 0; i < n; i++) free(inverse[i]);
            free(inverse);

            char cont;
            print_border(COLOR_BRIGHT_CYAN, 50);
            printf("%sBan co muon tiep tuc? (y/n): ", COLOR_BRIGHT_YELLOW);
            scanf(" %c", &cont);
            while (getchar() != '\n');
            print_border(COLOR_BRIGHT_CYAN, 50);
            if (cont == 'n' || cont == 'N') {
                for (int i = 0; i < n; i++) free(matrix[i]);
                free(matrix);
                system("cls");
                print_border(COLOR_BRIGHT_CYAN, 50);
                print_title(COLOR_BRIGHT_CYAN, "Cam on da su dung chuong trinh!", 50);
                print_border(COLOR_BRIGHT_CYAN, 50);
                return 0;
            }
            break;
        }
    }
    return 0;
}
