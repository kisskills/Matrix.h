#ifndef SRC_S21_MATRIX_H_
#define SRC_S21_MATRIX_H_

#define SUCCESS 1
#define FAILURE 0

#include <math.h>
#include <stdlib.h>
#include <float.h>

typedef struct matrix_struct {
    double** matrix;
    int rows;
    int columns;
} matrix_t;

typedef enum {
    OK = 0,
    INCORRECT_MATRIX = 1,
    CALC_ERROR = 2,
} res_code;

int s21_create_matrix(int rows, int columns, matrix_t *result);
void s21_remove_matrix(matrix_t *A);

int s21_eq_matrix(matrix_t *A, matrix_t *B);

int s21_sum_matrix(matrix_t *A, matrix_t *B, matrix_t *result);
int s21_sub_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_mult_number(matrix_t *A, double number, matrix_t *result);
int s21_mult_matrix(matrix_t *A, matrix_t *B, matrix_t *result);

int s21_transpose(matrix_t *A, matrix_t *result);

int s21_determinant(matrix_t *A, double *result);
int s21_calc_complements(matrix_t *A, matrix_t *result);

int s21_inverse_matrix(matrix_t *A, matrix_t *result);

void print_matrix(matrix_t A);
void add_matrix(matrix_t *A);
void init_matrix_by_number(double number, matrix_t *A);
int check_matrix(matrix_t *A);
int matrix_minor(matrix_t *A, int rows_number, int columns_number, matrix_t *result);
int check_zero_rows(matrix_t *A);
int check_zero_columns(matrix_t *A);
double determinant(matrix_t *A);
#endif  // SRC_S21_MATRIX_H_
