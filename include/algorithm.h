/**
 * \brief Function Prototype Library to Implement Sensor Fusion Algorithm
 * \author kaixia
 * \version 1.0
 * \date 2019-12-10
 */

#ifndef __algorithm_h__
#define __algorithm_h__


/**
 brief Structure for storing sensor
 details For the given some data of sensor which including number_of_sensor, value_of_each_sensor
 */
typedef struct sensor_list{
    int number_of_sensor;
    double value;
}sensor;

/**
 brief Structure for storing matrix
 details For the set of matrix data will store in the struct
 */
typedef struct
{
    int row;
    int column;
    double *data;
}matrix;

/**
 brief Function to initialize matrix
 details Initializing each new matrix when defined new matrix
 param[in] *matrix - it is a pointer to a struct of matrix
 param[in] row - it is a number of the row of the matrix
 param[in] colum - it is a number of the colum of the matrix
 return Success: true
 return Failure: false
 */
bool init_matrix(matrix *matrix, int row, int column);

/**
 brief Function to get the size of a matrix
 details getting a size of a matrix
 param[in] *matrix - it is a pointer to a struct of matrix
 return Success: size of the matrix
 return Failure: NULL
 */
int get_matrix_size(matrix *matrix);

/**
 brief Function to set matrix equal to zero
 details setting the every value in the matrix equal to zero
 param[in] *matrix - it is a pointer to a struct of matrix
 return Success: NULL
 return Failure: NULL
 */
void set_matrix_zeros(matrix *matrix);

/**
 brief Function to tell if the matrix is null
 details telling if the matrix is null
 param[in] *matrix - it is a pointer to a struct of matrix
 return Success: true
 return Failure: false
 */
bool is_null_matrix(matrix *matrix);

/**
 brief Function to destroy the matrix
 details destroying matrix that will no longer use
 param[in] *matrix - it is a pointer to a struct of matrix
 return Success: NULL
 return Failure: NULL
 */
void destroy_matrix(matrix *matrix);

/**
 brief Function to compute the second matrix norm
 details computing the second norm of the matrix
 param[in] *matrix - it is a pointer to a struct of matrix
 return Success: the value of the matrix norm
 return Failure: NULL
 */
double matrix_norm2(matrix *matrix);

/**
 brief Function to copy the matrix
 details copying the value of matrix_a to the value of matrix_b
 param[in] *matrix_a - it is a pointer to a struct of matrix a
 param[in] *matrix_b - it is a pointer to a struct of matrix b
 return Success: NULL
 return Failure: NULL
 */
void copy_matrix(matrix *matrix_a, matrix *matrix_b);

/**
 brief Function to decompose the matrix
 details processing the matrix using QR-decomposition method, matrix_a is the original matirx
 param[in] *a - it is a pointer to a struct of matrix a
 param[in] *q - it is a pointer to a struct of matrix q
 param[in] *r - it is a pointer to a struct of matrix r
 return Success: NULL
 return Failure: NULL
 */
void qr(matrix *a, matrix *q, matrix *r);

/**
 brief Function to sort the eigenvalues
 details sorting all eigenvalues, when falg is equal to 0 that means desc, it is equal to 1 that means asce
 param[in] *eigen_value - it is a pointer to a struct of matrix eigen_value
 param[in] flag - it is a flag to change the sorted order
 return Success: ture
 return Failure: false
 */
bool sort_eigen_values(matrix *eigen_value, int flag);

/**
 brief Function to multiply two matrixs
 details multiplying the values of two different matrixs a,b and getting the result matrix c
 param[in] *c - it is a pointer to a struct of matrix c
 param[in] *a - it is a pointer to a struct of matrix a
 param[in] *b - it is a pointer to a struct of matrix b
 return Success: ture
 return Failure: false
 */
bool matrix_mul_matrix(matrix *c, matrix *a, matrix *b);

/**
 brief Function to comput the eigenvectors
 details calculating all eigenvetors through all eigenvalues we got before, and storing in matrix a
 param[in] *a - it is a pointer to a struct of matrix a
 param[in] *eigen_vector - it is a pointer to a struct of matrix eigen_vector
 param[in] *eigen_values - it is a pointer to a struct of matrix eigen_values
 return Success: NULL
 return Failure: NULL
 */
void calculate_eigenvector(matrix *a, matrix *eigen_vector, matrix *eigen_values);

/**
 brief Function to get support degree matrix
 details calculating the support degree of every two sensors and getting the support degree matrix
 param[in] sensor - it is a struct of storing all sensor data
 param[in] n - it is a number of all sensors
 return Success: a pointer to a two dimentional array of support degree matrix
 return Failure: NULL
 */
double **calculate_sd_matrix(struct sensor_list sensor[],int n);

/**
 brief Function to transpose eigenvectors
 details generating all transposed eigenvectors
 param[in] *matrix - it is a pointer to a struct of matrix matrix
 return Success: a pointer to a two dimentional array of d
 return Failure: NULL
 */
double **generate_eigenvector_transpose(matrix *matrix);

/**
 brief Function to calculate the value of alpha_k
 details generating an array of all values of contribution rate of kth principal component
 param[in] *matrix - it is a pointer to a struct of matrix matrix
 param[in] n - it is a number of all sensors
 return Success: a pointer to an array of aplha_k
 return Failure: NULL
 */
double *calculate_contribution_rate_of_kth_principal_component(matrix *matrix,int n);

/**
 brief Function to calculate the value of z
 details computing an array of all values of integrated support degree score
 param[in] *matrix - it is a pointer to a struct of matrix matrix
 param[in] *alpha_k - it is a pointer to an array of aplha_k
 param[in] n - it is a number of all sensors
 param[in] m - it is a value that means we will select the first m principal components later
 return Success: a pointer to an array of integrated support degree score
 return Failure: NULL
 */
double *compute_integrated_support_degree_score(matrix *y,double *alpha_k,int n,int m);

/**
 brief Function to disregard wrong data
 details disregarding invalid data, setting the invalid sensor equal to 0.0 and getting an array of new z after disregarding
 param[in] *z - it is a pointer to an array of all values of integrated support degree score
 param[in] n - it is a number of all sensors
 param[in] q - it is a changeable value that can limit the valid data
 return Success: a pointer to an array of new integrated support degree score after disgarding wrong sensor data
 return Failure: NULL
 */
double *disregard_wrong_data(double *z,int n,double q);

/**
 brief Function to compute weight coefficient
 details computing the weight coefficient of every valid sensor and getting an array of weight
 param[in] *z_after - it is a pointer to an array new integrated support degree score after disgarding wrong sensor data
 param[in] n - it is a number of all sensors
 return Success: a pointer to an array of weight coefficient of every valid sensor
 return Failure: NULL
 */
double *compute_weight_coefficient(double *z_after,int n);

/**
 brief Function to compute fused output
 details getting the value of fused output as the fianl result
 param[in] sensor - it is a struct of storing all sensor data
 param[in] *weight - it is a pointer to an array of weight coefficient of every valid sensor
 param[in] n - it is a number of all sensors
 return Success: the value of fused output as the fianl result
 return Failure: NULL
 */
double compute_fused_output(struct sensor_list sensor[],double *weight,int n);


#endif //__algorithm_h__

