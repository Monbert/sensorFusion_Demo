/**
 * \brief Function Definition Library to Implement Data_input and Data_output
 * \details this file implements the input and output sensor data and 
 *    puts data into the sensorfuion algorithm
 * \author kaixia
 * \version 1.0
 * \date 2019-12-13
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "../include/algorithm.h"
#include "../include/data.h"

/**
 * \brief Function to output data into an .xls file
 * \details write the final outputs fused_value to the output file
 * \param[in] output_1 - the final fused_value
 * \param[in] output_2 - the final fused_value
 * \return Success: NULL
 * \return Failure: NULL
 */
void output_data(double output_1,double output_2){
    FILE *fp;
    fp=fopen("./data/output/outputdata.xls","w");
    fseek(fp,0,SEEK_SET);
    fprintf(fp,"time_24h\tfused_value\n");
    fprintf(fp,"%3.1f\t",13.2);
    fprintf(fp,"%f\n",output_1);
    fprintf(fp,"%3.1f\t",13.5);
    fprintf(fp,"%f\n",output_2);
}

/**
 * \brief Function to read the data from the input .csv file
 * \details read the initial data from sensors and store them into a strcut variable
 * \parameter[in] - NULL
 * \return Success: the initial address of the struct variable sensor
 * \return Failure: NULL
 */
struct sensor_list *input_data(){
    FILE *fp;
    int i, j;
    double D[16][3] = {};
    fp=fopen("./data/input/inputdata.csv","r");
    if (fp == NULL){
        printf("error! There is no inputdata.csv file, please "
            "check if it is existed.\n");
        exit(0);
    }
    fseek(fp,25L,SEEK_SET);
    for(i=0; i<16; i++)
        for(j= 0 ; j < 3; j++){
            fscanf(fp,"%lf",&D[i][j]);
            fseek(fp, 1L, SEEK_CUR);
        }
    static struct sensor_list sensor[16];
    for(i=0; i<16; i++)
        for(j=0; j<3; j++){
            sensor[i].number_of_sensor=i;
            sensor[i].value=D[i][2];
        }
    fclose(fp);
    return sensor;
}

/**
 * \brief Function to set input sensor data into algorithm flow
 * \details repeatign to set input sensor data into algorithm flow and get results
 * \parameter[in] sensor - it is a struct of storing all sensor data
 * \return Success: the value of output sensor data
 * \return Failure: NULL
 */
double repeat(struct sensor_list sensor[]){
    int i,j,k,m,n,row,col,num;
    n = 8;
    m = 0;
    row = n;
    col = n;
    num = row * col;
    double sum_varphi_m = 0.0;
    double p = 0.85;
    double q = 0.7;
    double ** d;
    double ** eigenvectors_transpose_array;
    double * alpha_k;
    double * z;
    double * z_after;
    double * weight;
    double output;
    const unsigned number_of_iterations = 80;
    matrix mymatrix, temp, temp_q, temp_r,eigen_value, eigen_vectors, y;
    
    init_matrix(&mymatrix, row, col);
    init_matrix(&temp, row, col);
    init_matrix(&y, row, col);
    set_matrix_zeros(&temp);
    set_matrix_zeros(&mymatrix);
    set_matrix_zeros(&y);
    
    d = calculate_sd_matrix(sensor, n);
    for(i=0, k=0; i<n&&k<num; i++){
        for(j=0; j<n; j++){
            mymatrix.data[k] = d[i][j];
            k++;
        }
    }
    copy_matrix(&mymatrix, &temp);
    init_matrix(&temp_q, mymatrix.row, mymatrix.column);
    init_matrix(&temp_r, mymatrix.row, mymatrix.column);
    init_matrix(&eigen_value, mymatrix.row, 1);

    for (k = 0; k < number_of_iterations; ++k){
        qr_decomposition(&temp, &temp_q, &temp_r);
        matrix_mul_matrix(&temp, &temp_r, &temp_q);
    }
    for (k = 0; k < temp.column; ++k){
        eigen_value.data[k] = temp.data[k * temp.column + k];
    }
    sort_eigen_values(&eigen_value, 0);
    
    calculate_eigenvector(&mymatrix, &temp_q, &eigen_value);
    eigenvectors_transpose_array = generate_eigenvector_transpose(&temp_q);
    init_matrix(&eigen_vectors, row, col);
    for(i=0, k=0; i<n&&k<num; i++){
        for(j=0; j<n; j++){
            eigen_vectors.data[k] = eigenvectors_transpose_array[i][j];
            k++;
        }
    }
    matrix_mul_matrix(&y, &eigen_vectors, &mymatrix);
    alpha_k = calculate_contribution_rate_of_kth_principal_component(&eigen_value, n);
    for(i=0; i<n; i++){
        sum_varphi_m += alpha_k[i];
        m = i + 1;
        if(sum_varphi_m > p) break;
    }
    z = compute_integrated_support_degree_score(&y,alpha_k, n, m);
    z_after = disregard_wrong_data(z, n, q);
    weight = compute_weight_coefficient(z_after, n);
    output = compute_fused_output(sensor,weight, n);
    
    destroy_matrix(&eigen_value);
    destroy_matrix(&mymatrix);
    destroy_matrix(&temp);
    destroy_matrix(&temp_q);
    destroy_matrix(&temp_r);
    destroy_matrix(&y);
    for (i = 0; i < n; ++i){
        free(d[i]);
    }
    free(d);
    for (i = 0; i < n; ++i){
        free(eigenvectors_transpose_array[i]);
    }
    free(eigenvectors_transpose_array);
    free(alpha_k);
    free(z);
    free(z_after);
    free(weight);
    return output;
}
