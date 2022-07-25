#include "s21_matrix.h"

void init_matrix_by_number(double number, matrix_t *A) {
    for (int i = 0; i < A->rows; ++i) {
        for (int j = 0; j < A->columns; ++j) {
            A->matrix[i][j] = number;
            number += 1.0;
        }
    }
}

int check_matrix(matrix_t *A) {
    int res = OK;

    if (!A || A->rows < 1 || A->rows < 1) {
        res = INCORRECT_MATRIX;
    }

    return res;
}

int matrix_minor(matrix_t *A, int rows_number, int columns_number, matrix_t *result) {
    if (check_matrix(A)) return INCORRECT_MATRIX;
    if (A->rows != A->columns) return CALC_ERROR;

    int res = s21_create_matrix(A->rows - 1, A->columns - 1, result);
    if (res == OK) {
        for (int i = 0, i_old = 0; i < result->rows; i++, i_old++) {
            if (i_old == rows_number) i_old++;
            for (int j = 0, j_old = 0; j < result->columns; j++, j_old++) {
                if (j_old == columns_number) j_old++;
                result->matrix[i][j] = A->matrix[i_old][j_old];
            }
        }
    }
    return res;
}

int check_zero_rows(matrix_t *A) {
    int res = 0;

    for (int i = 0; i < A->rows; i++) {
        for (int j = 0; j < A->columns; j++) {
            if (A->matrix[i][j] != 0.0) break;
            if ((j == A->columns - 1) && (A->matrix[i][j] == 0)) res = 1;
        }
        if (res) break;
    }
    return res;
}

int check_zero_columns(matrix_t *A) {
    int res = 0;

    for (int j = 0; j < A->columns; j++) {
        for (int i = 0; i < A->rows; i++) {
            if (A->matrix[i][j] != 0.0) break;
            if ((i == A->rows - 1) && (A->matrix[i][j] == 0)) res = 1;
        }
        if (res) break;
    }
    return res;
}

double determinant(matrix_t *A) {
    if (A->rows == 1) return A->matrix[0][0];
    if (A->rows == 2) return A->matrix[0][0] * A->matrix[1][1] - A->matrix[0][1] * A->matrix[1][0];

    double res = 0.0;
    matrix_t minor_matrix = {0};
    for (int k = 0; k < A->columns; k++) {
        matrix_minor(A, 0, k, &minor_matrix);
        if ((k + 2) % 2)
            res -= (A->matrix[0][k] * determinant(&minor_matrix));
        else
            res += (A->matrix[0][k] * determinant(&minor_matrix));
        s21_remove_matrix(&minor_matrix);
    }
    return res;
}

/* functions for testing */
// void print_matrix(matrix_t A) {
//     printf("rows = %d\tcolumns = %d\n", A.rows, A.columns);
//     for (int i = 0; i < A.rows; i++) {
//         for (int j = 0; j < A.columns; j++) {
//             printf("%7.3lf ", A.matrix[i][j]);
//         }
//         printf("\b\n");
//     }
// }

// void add_matrix(matrix_t *A) {
//     for (int i = 0; i < A->rows; i++) {
//         for (int j = 0; j < A->columns; j++) {
//             scanf("%lf", &(A->matrix[i][j]));
//         }
//     }
// }
