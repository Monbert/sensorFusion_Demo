/**
 * \brief Function Definition Library to Implement Sensor Fusion Algorithm
 * \details this file implements the sensorfusion core computing steps, 
 *    which we call sensorfusion algorithm
 * \author kaixia
 * \version 1.0
 * \date 2019-12-10
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "../include/algorithm.h"

/**
 * \brief Function to initialize matrix
 * \details Initializing each new matrix when defined new matrix
 * \param[in] *matrix - it is a pointer to a struct of matrix
 * \param[in] row - it is a number of the row of the matrix
 * \param[in] column - it is a number of the column of the matrix
 * \return Success: true
 * \return Failure: false
 */
bool init_matrix(matrix *matrix, int row, int column){
    int matrix_size = row*column*sizeof(double);
    if (matrix_size <= 0){
        return false;
    }
    matrix->data = (double*)malloc(matrix_size);
    if (matrix->data){
        matrix->row = row;
        matrix->column = column;
        return true;
    }
    else{
        matrix->row = 0;
        matrix->column = 0;
        return false;
    }
}

/**
 * \brief Function to get the size of a matrix
 * \details getting a size of a matrix
 * \param[in] *matrix - it is a pointer to a struct of matrix
 * \return Success: size of the matrix
 * \return Failure: NULL
 */
int get_matrix_size(matrix *matrix){
    return matrix->row*matrix->column;
}

/**
 * \brief Function to set matrix equal to zero
 * \details setting the every value in the matrix equal to zero
 * \param[in] *matrix - it is a pointer to a struct of matrix
 * \return Success: NULL
 * \return Failure: NULL
 */
void set_matrix_zeros(matrix *matrix){
    int matrix_num = get_matrix_size(matrix);
    for (int i = 0; i < matrix_num; i++){
        matrix->data[i] = 0;
    }
}

/**
 * \brief Function to tell if the matrix is null
 * \details telling if the matrix is null
 * \param[in] *matrix - it is a pointer to a struct of matrix
 * \return Success: true
 * \return Failure: false
 */
bool is_null_matrix(matrix *matrix){
    int matrix_num =get_matrix_size(matrix);
    if ((matrix_num <= 0) || (matrix->data == NULL)){
        return true;
    }
    else{
        return false;
    }
}

/**
 * \brief Function to destroy the matrix
 * \details destroying matrix that will no longer use
 * \param[in] *matrix - it is a pointer to a struct of matrix
 * \return Success: NULL
 * \return Failure: NULL
 */
void destroy_matrix(matrix *matrix){
    if (!is_null_matrix(matrix)){
        matrix->data = NULL;
        matrix->row = 0;
        matrix->column = 0;
        free(matrix->data);
    }
}

/**
 * \brief Function to compute the matrix norm
 * \details computing the matrix norm and return this 
 *    value to the function which needed
 * \param[in] *matrix - it is a pointer to a struct of matrix
 * \return Success: the value of the matrix norm
 * \return Failure: NULL
 */
double get_matrix_norm(matrix *matrix){
    double matrix_norm = 0.0;
    int matrix_num = get_matrix_size(matrix);
    for (int i = 0; i < matrix_num; i++){
        matrix_norm+=(matrix->data[i]) * (matrix->data[i]);
    }
    matrix_norm = (double)sqrt(matrix_norm);
    return matrix_norm;
}

/**
 * \brief Function to copy the matrix
 * \details copying the value of matrix_a to the value of matrix_b
 * \param[in] *matrix_a - it is a pointer to a struct of matrix a
 * \param[in] *matrix_b - it is a pointer to a struct of matrix b
 * \return Success: NULL
 * \return Failure: NULL
 */
void copy_matrix(matrix *matrix_a, matrix *matrix_b){
    if (matrix_b->row != matrix_a->row){
        matrix_b->row = matrix_a->row;
    }
    if (matrix_b->column != matrix_a->column){
        matrix_b->column = matrix_a->column;
    }
    int size_a = get_matrix_size(matrix_a);
    for (int i = 0; i < size_a; i++){
        matrix_b->data[i] = matrix_a->data[i];
    }
}

/**
 * \brief Function to go through QR decomposition
 * \details processing the matrix using QR decomposition method, 
 *    matrix_a is the original matirx
 * \param[in] *a - it is a pointer to a struct of matrix a, 
 *    which is a original matrix
 * \param[in] *q - it is a pointer to a struct of matrix q, 
 *    which is an orthogonal matrix
 * \param[in] *r - it is a pointer to a struct of matrix r, 
 *    which is an upper triangular matrix
 * \return Success: NULL
 * \return Failure: NULL
 */
void qr_decomposition(matrix *a, matrix *q, matrix *r){
    matrix col_a, col_q;
    init_matrix(&col_a, a->row, 1);
    set_matrix_zeros(&col_a);
    init_matrix(&col_q, a->row, 1);
    set_matrix_zeros(&col_q);
    if (a->row != a->column){
        printf("a is not a square matrix!");
    }
    int a_size = get_matrix_size(a);
    int q_size = get_matrix_size(q);
    int r_size = get_matrix_size(r);
    if (q_size != a_size){
        destroy_matrix(q);
        init_matrix(q, a->row, a->column);
        set_matrix_zeros(q);
    }
    else{
        q->row = a->row;
        q->column = a->column;
        set_matrix_zeros(q);
    }
    if (r_size != a_size){
        destroy_matrix(r);
        init_matrix(r, a->row, a->column);
        set_matrix_zeros(r);
    }
    else{
        r->row = a->row;
        r->column = r->column;
        set_matrix_zeros(r);
    }
    for (int j = 0; j < a->column; j++){
        for (int i = 0; i < a->column; i++){
            col_a.data[i] = a->data[i * a->column + j];
            col_q.data[i] = a->data[i * a->column + j];
        }
        for (int k = 0; k < j; k++){
            r->data[k * r->column + j] = 0;
            for(int i1 = 0; i1 < col_a.row; i1++){
                r->data[k * r->column + j] += col_a.data[i1] * 
                q->data[i1 * q->column + k];
            }
            for (int i2 = 0; i2 < a->column; i2++){
                col_q.data[i2] -= r->data[k * r->column + j] * 
                q->data[i2 * q->column + k];
            }
        }
        double temp = get_matrix_norm(&col_q);
        r->data[j * r->column + j] = temp;
        for (int i3 = 0; i3 < q->column; i3++){
            q->data[i3 * q->column + j] = col_q.data[i3] / temp;
        }
    }
    destroy_matrix(&col_a);
    destroy_matrix(&col_q);
}

/**
 * \brief Function to sort the eigenvalues
 * \details sorting all eigenvalues, when flag is equal to 0 that means desc, 
 *    it is equal to 1 that means asce
 * \param[in] *eigen_value - it is a pointer to a struct of matrix eigen_value
 * \param[in] flag - it is a flag to change the sorted order
 * \return Success: ture
 * \return Failure: false
 */
bool sort_eigen_values(matrix *eigen_value, int flag){
    int size = get_matrix_size(eigen_value);
    for (int i = 0; i < size - 1; i++){
        int k = i;
        for (int j = i + 1; j < size; j++){
            if (flag == 1){
                if (eigen_value->data[k] > eigen_value->data[j]){
                    k = j;
                }
            }
            else if (flag == 0){
                if (eigen_value->data[k] < eigen_value->data[j]){
                    k = j;
                }
            }
            else{
                return false;
            }
        }
        if (k != i){
            double temp;
            temp = eigen_value->data[i];
            eigen_value->data[i] = eigen_value->data[k];
            eigen_value->data[k] = temp;
        }
    }
    return true;
}

/**
 * \brief Function to multiply two matrixs
 * \details multiplying the values of two different matrixs a,b and 
 *    getting the result matrix c
 * \param[in] *c - it is a pointer to a struct of matrix c
 * \param[in] *a - it is a pointer to a struct of matrix a
 * \param[in] *b - it is a pointer to a struct of matrix b
 * \return Success: ture
 * \return Failure: false
 */
bool matrix_mul_matrix(matrix *c, matrix *a, matrix *b){
    if ((is_null_matrix(a)) || (is_null_matrix(b))){
        return false;
    }
    int a_col = a->column;
    int b_row = b->row;
    init_matrix(c, a->row, b->column);
    set_matrix_zeros(c);
    if (a_col != b_row){
        printf("a_col!=b_row!");
        return false;
    }
    for (int i = 0; i < a->row; i++){
        for (int j = 0; j < b->column; j++)
        {
            for (int k = 0; k < a->column; k++){
                c->data[i*c->row + j] += a->data[i*a->row + k] * 
                b->data[k*b->column + j];
            }
        }
    }
    return true;
}

/**
 * \brief Function to comput the eigenvectors
 * \details calculating all eigenvetors through all eigenvalues we got before, 
 *    and storing in matrix a
 * \param[in] *a - it is a pointer to a struct of matrix a
 * \param[in] *eigen_vector - it is a pointer to a struct of matrix eigen_vector
 * \param[in] *eigen_values - it is a pointer to a struct of matrix eigen_values
 * \return Success: NULL
 * \return Failure: NULL
 */
void calculate_eigenvector(matrix *a, matrix *eigen_vector, matrix *eigen_values){
    int num = a->column;
    double eigen_value;
    matrix temp;
    init_matrix(&temp, a->row, a->column);
    for (int count = 0; count < num; count++){
        eigen_value = eigen_values->data[count];
        copy_matrix(a, &temp);
        for (int i = 0; i < temp.row; i++){
            temp.data[i * temp.column + i] -= eigen_value;
        }
        for (int i = 0; i < temp.row - 1; i++){
            double coe = temp.data[i * temp.column + i];
            for (int j = i; j < temp.column; j++){
                temp.data[i * temp.column + j] /= coe;
            }
            for (int i1 = i + 1; i1 < temp.row; i1++){
                coe = temp.data[i1 * temp.column + i];
                for (int j1 = i; j1 < temp.column; j1++){
                    temp.data[i1 * temp.column + j1] -= coe * 
                    temp.data[i * temp.column + j1];
                }
            }
        }
        double sum1 = eigen_vector->data[(eigen_vector->row - 1) * 
            eigen_vector->column + count] = 1;
        for (int i2 = temp.row - 2; i2 >= 0; i2--){
            double sum2 = 0;
            for (int j2 = i2 + 1; j2 < temp.column; j2++){
                sum2 += temp.data[i2 * temp.column + j2] * 
                eigen_vector->data[j2 * eigen_vector->column + count];
            }
            sum2 = -sum2 / temp.data[i2 * temp.column + i2];
            sum1 += sum2 * sum2;
            eigen_vector->data[i2 * eigen_vector->column + count] = sum2;
        }
        sum1 = sqrt(sum1);
        for (int i = 0; i < eigen_vector->row; i++)
        {
            eigen_vector->data[i * eigen_vector->column + count] /= sum1;
        }
    }
    destroy_matrix(&temp);
}

/**
 * \brief Function to get support degree matrix
 * \details calculating the support degree of every two sensors and 
 *    getting the support degree matrix
 * \param[in] sensor - it is a struct of storing all sensor data
 * \param[in] n - it is a number of all sensors
 * \return Success: a pointer to a two dimentional array of support degree matrix
 * \return Failure: NULL
 */
double **calculate_sd_matrix(struct sensor_list sensor[],int n){
    int i, j;
    double **d = (double **)malloc(sizeof(double *) * n);
    for (i = 0; i < n; ++i){
        d[i] = (double *)malloc(sizeof(double) * n);
    }
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            double temp = sensor[i].value - sensor[j].value;
            if(temp<0){
                temp = -temp;
                d[i][j] = exp(-(temp));
            }
            else{
                d[i][j] = exp(-(temp));
            }
        }
    }
    return d;
}

/**
 * \brief Function to transpose eigenvectors
 * \details generating all transposed eigenvectors
 * \param[in] *matrix - it is a pointer to a struct of matrix matrix
 * \return Success: a pointer to a two dimentional array of d
 * \return Failure: NULL
 */
double **generate_eigenvector_transpose(matrix *matrix){
    int i,j;
    int n = matrix->row;
    double **d = (double **)malloc(sizeof(double *) * n);
    for (i = 0; i < n; ++i){
        d[i] = (double *)malloc(sizeof(double) * n);
    }
    
    for(i=0; i<n; i++){
        for(j=0; j<n; j++){
            d[i][j] = matrix->data[i+n*j];
        }
    }
    return d;
}

/**
 * \brief Function to calculate the value of alpha_k
 * \details generating an array of all values of contribution rate of 
 *    kth principal component
 * \param[in] *matrix - it is a pointer to a struct of matrix matrix
 * \param[in] n - it is a number of all sensors
 * \return Success: a pointer to an array of aplha_k
 * \return Failure: NULL
 */
double *calculate_contribution_rate_of_kth_principal_component(
    matrix *matrix, int n){
    double *alpha_k;
    double sum = 0.0;
    alpha_k = (double*)malloc(sizeof(double)*n);
    for(int i=0; i<n; i++){
        sum += matrix->data[i];
    }
    for(int i=0; i<n; i++){
        alpha_k[i] = matrix->data[i] / sum;
    }
    return alpha_k;
}

/**
 * \brief Function to calculate the value of z
 * \details computing an array of all values of integrated support degree score
 * \param[in] *matrix - it is a pointer to a struct of matrix matrix
 * \param[in] *alpha_k - it is a pointer to an array of aplha_k
 * \param[in] n - it is a number of all sensors
 * \param[in] m - it is a value that means we will select the first 
 *    m principal components later
 * \return Success: a pointer to an array of integrated support degree score
 * \return Failure: NULL
 */
double *compute_integrated_support_degree_score(
    matrix *y,double *alpha_k,int n,int m){
    double *z;
    int i,j;
    z = (double*)malloc(sizeof(double)*n);
    for(i=0; i<n; i++){
        for(j=0; j<m; j++){
            z[i] += y->data[n*j+i] * alpha_k[j];
        }
    }
    return z;
}

/**
 * \brief Function to disregard wrong data
 * \details disregarding invalid data, setting the invalid sensor equal to 0.0 and 
 *    getting an array of new z after disregarding
 * \param[in] *z - it is a pointer to an array of all values of integrated support degree score
 * \param[in] n - it is a number of all sensors
 * \param[in] q - it is a changeable value that can limit the valid data
 * \return Success: a pointer to an array of new integrated support degree score 
 *    after disgarding wrong sensor data
 * \return Failure: NULL
 */
double *disregard_wrong_data(double *z,int n,double q){
    int i;
    double sum = 0.0;
    double *z_after;
    z_after = (double*)malloc(sizeof(double)*n);
    
    for(i=0; i<n; i++){
        sum += z[i];
    }
    for(i=0; i<n; i++){
        if((z[i]*z[i]) < (sum*sum/n/n*q*q)){
            z_after[i] = 0.0;
            printf("The value of sensor number[%d] is outside a prescribed range\n",i+1);
        }
        else{
            z_after[i] = z[i];
        }
    }
    return z_after;
}

/**
 * \brief Function to compute weight coefficient
 * \details computing the weight coefficient of every valid sensor and 
 *    getting an array of weight
 * \param[in] *z_after - it is a pointer to an array 
 *    new integrated support degree score after disgarding wrong sensor data
 * \param[in] n - it is a number of all sensors
 * \return Success: a pointer to an array of weight coefficient of every valid sensor
 * \return Failure: NULL
 */
double *compute_weight_coefficient(double *z_after,int n){
    double *weight;
    double sum_z_after =0.0;
    double sum_weight =0.0;
    int i;
    weight = (double*)malloc(sizeof(double)*n);
    for(i=0; i<n; i++){
        sum_z_after += z_after[i];
    }
    for(i=0; i<n; i++){
        weight[i] = z_after[i] / sum_z_after;
        sum_weight += weight[i];
    }
    return weight;
}

/**
 * \brief Function to compute fused output
 * \details getting the value of fused output as the fianl result
 * \param[in] sensor - it is a struct of storing all sensor data
 * \param[in] *weight - it is a pointer to an array of 
 *    weight coefficient of every valid sensor
 * \param[in] n - it is a number of all sensors
 * \return Success: the value of fused output as the fianl result
 * \return Failure: NULL
 */
double compute_fused_output(struct sensor_list sensor[],double *weight,int n){
    double x_output=0.0;
    int i;
    for(i=0; i<n; i++){
        x_output += weight[i]*sensor[i].value;
    }
    return x_output;
}

