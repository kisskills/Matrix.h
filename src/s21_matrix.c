#include "s21_matrix.h"

int s21_create_matrix(int rows, int columns, matrix_t *result) {
    int res = OK;

    if ((rows < 1) || (columns < 1)) return INCORRECT_MATRIX;

    result->matrix = calloc(1, rows * columns * sizeof(double) + rows * sizeof(double*));
    if (result->matrix) {
        double *matrix_ptr = (double*) (result->matrix + rows);
        for (int i = 0; i < rows; i++)
            result->matrix[i] = matrix_ptr + i * columns;
    } else {
        res = CALC_ERROR;
    }

    result->rows = rows;
    result->columns = columns;

    return res;
}

void s21_remove_matrix(matrix_t *A) {
    if (A->matrix)
        free(A->matrix);

    A->rows = 0;
    A->columns = 0;
    A->matrix = NULL;
}

int s21_eq_matrix(matrix_t *A, matrix_t *B) {
    int res = SUCCESS;

    if ((A->rows == B->rows) && (A->columns == B->columns)) {
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                if ((A->matrix[i][j] < 0 && B->matrix[i][j] > 0) || \
                (A->matrix[i][j] > 0 && B->matrix[i][j] < 0)) res = FAILURE;
                else if (fabs(A->matrix[i][j] - B->matrix[i][j]) > 1e-7) res = FAILURE;
                if (!res) break;
            }
            if (!res) break;
        }
    } else {
        res = FAILURE;
    }
    return res;
}

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int res = OK;

    if ((A->rows == B->rows) && (A->columns == B->columns)) {
        res = s21_create_matrix(A->rows, A->columns, result);
        if (res == OK) {
            for (int i = 0; i < result->rows; i++)
                for (int j = 0; j < result->columns; j++)
                    result->matrix[i][j] = A->matrix[i][j] + B->matrix[i][j];
        }
    } else {
        res = CALC_ERROR;
    }

    return res;
}

int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    int res = OK;

    if (A->rows == B->rows && A->columns == B->columns) {
        res = s21_create_matrix(A->rows, A->columns, result);
        if (res == OK) {
            for (int i = 0; i < A->rows; i++)
                for (int j = 0; j < A->columns; j++)
                    result->matrix[i][j] = A->matrix[i][j] - B->matrix[i][j];
        }
    } else {
        res = CALC_ERROR;
    }

    return res;
}

int s21_mult_number(matrix_t *A, double number, matrix_t *result) {
    int res = s21_create_matrix(A->rows, A->columns, result);
    if (res == OK) {
        for (int i = 0; i < A->rows; i++)
            for (int j = 0; j < A->columns; j++)
                result->matrix[i][j] = A->matrix[i][j] * number;
    }

    return res;
}

int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result) {
    if (check_matrix(A) || check_matrix(B)) return INCORRECT_MATRIX;
    if (A->columns != B->rows) return CALC_ERROR;
    int res = s21_create_matrix(A->rows, B->columns, result);

    if (res == OK) {
        for (int i = 0; i < A->rows; i++)
            for (int k = 0; k < A->columns; k++)
                for (int j = 0; j < B->columns; j++)
                    result->matrix[i][j] += A->matrix[i][k] * B->matrix[k][j];
    }

    return OK;
}

int s21_transpose(matrix_t *A, matrix_t *result) {
    if (check_matrix(A)) return INCORRECT_MATRIX;
    int res = s21_create_matrix(A->columns, A->rows, result);

    if (res == OK) {
        for (int i = 0; i < A->columns; i++)
            for (int j = 0; j < A->rows; j++)
                result->matrix[i][j] = A->matrix[j][i];
    }

    return res;
}

int s21_calc_complements(matrix_t *A, matrix_t *result) {
    if (check_matrix(A)) return INCORRECT_MATRIX;
    if (A->rows != A->columns || A->rows == 1) return CALC_ERROR;

    double det = 0.0;
    int res = s21_create_matrix(A->rows, A->columns, result);
    matrix_t minor_matrix = {0};

    if (res == OK) {
        for (int i = 0; i < A->rows; i++) {
            for (int j = 0; j < A->columns; j++) {
                res = matrix_minor(A, i, j, &minor_matrix);
                if (res == OK) res = s21_determinant(&minor_matrix, &det);
                if (res == OK) {
                    if ((i + j) % 2)
                        result->matrix[i][j] -= det;
                    else
                        result->matrix[i][j] += det;
                }
		s21_remove_matrix(&minor_matrix);
                if (res != OK) break;
            }
            if (res != OK) break;
        }
    }

    return res;
}

int s21_determinant(matrix_t *A, double *result) {
    if (check_matrix(A)) return INCORRECT_MATRIX;
    if (A->rows != A->columns) return CALC_ERROR;

    if (check_zero_rows(A) || check_zero_columns(A)) {
        *result = 0.0;
    } else {
        *result = determinant(A);
    }

    return OK;
}

int s21_inverse_matrix(matrix_t *A, matrix_t *result) {
    int res = OK;
    double det = 0;

    res = s21_determinant(A, &det);
    if (res != OK) return res;
    if (A->rows == 1) {
        s21_create_matrix(A->rows, A->columns, result);
        result->matrix[0][0] = 1.0 / A->matrix[0][0];
    } else if (fabs(det) > 1e-6) {
	matrix_t calc_matrix = {0};
	matrix_t trans_matrix = {0};
        s21_create_matrix(A->rows, A->columns, &calc_matrix);
        s21_calc_complements(A, &calc_matrix);
        s21_transpose(&calc_matrix, &trans_matrix);
	s21_remove_matrix(&calc_matrix);

        s21_mult_number(&trans_matrix, 1.0 / det, result);

        s21_remove_matrix(&trans_matrix);
    } else {
        res = CALC_ERROR;
    }

    return res;
}
