/**
 * \brief test of each individual function of our program
 * \details this file implements the sensorfusion core computing steps,
 *    which we call sensorfusion algorithm, and compare each output with our
 *    expected output to see whether individual component fullfil its job
 * \author kaixia
 * \version 1.0
 * \date 2019-12-17
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "../../include/algorithm.h"
#include "../../include/data.h"

/**
 * \brief  to test the correctness of our whole algorithm and 
 *    each individual function component, and do some initialization and finishing touches
 * \details to manually compute each expected outputs of each function and compare it with 
 *    our actual outputs from the function and set some evaluation metrics to see 
 *    whether each of our function meets our standard
 * \param[in]  - the initial address of the struct variable sensor_list sensor_1
 * \return Success: true: each of function output meets our expected output and return a convincing result
 * \return Failure: some of the function output don't meet our expected output or 
 *    the process don't run normally ,or the expected final output is not what we get from the main program.
 */
double test_sensordataset1(struct sensor_list sensor[]){
    int i,j,k,m,n,row,col,num;
    n = 8;
    m = 0;
    row = n;
    col = n;
    num = row*col;
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
    matrix mymatrix,temp,temp_q,temp_r, eigen_value, eigen_vectors, y;
    init_matrix(&mymatrix, row, col);

    printf("**testing init_matrix function**\n");
    printf("**param[in]- a pointer to a matrix**\n");
    printf("**param[in]- the row and column of the matrix**\n");
    printf("**expected result:if it run successfully:row=8,column=8**\n");
    printf("**expected result:if it fails:row or column not equal to 8**\n");
    printf("**factual result:row=%d\tcolumn=%d**\n",mymatrix.row,mymatrix.column);
    if(mymatrix.column==8&&mymatrix.row==8){
        printf("**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("**test result of this function:FAIL**\n\n\n");
    }
    int size=get_matrix_size(&mymatrix);

    printf("**testing get_matrix_size function**\n");
    printf("**param[in]- a pointer to a matrix**\n");
    printf("**param[out]-size of the matrix**\n");
    printf("**expected result:if it run successfully:size=64**\n");
    printf("**expected result:if it fails:size not equal to 64**\n");
    printf("**Factual result: %d**\n",size);
    if(size==64){
        printf("**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("**test result of this function:FAIL**\n\n\n");
    }
    init_matrix(&temp, row, col);
    init_matrix(&y, row, col);
    set_matrix_zeros(&temp);
    set_matrix_zeros(&mymatrix);
    set_matrix_zeros(&y);

    int mymatrix_number=0;
    printf("**testing set_matrix zeros function**\n");
    printf("**param[in]- a pointer to a matrix**\n");
    printf("**param[out]-NULL**\n");
    printf("**expected result:if it run successfully:mymatrix.data[i]==0**\n");
    printf("**expected result:if it fails:have mymatrix.data[i]!=0**\n");
    printf("**Factual result:\n");
    for(i=0;i<8;i++){
        printf("%7.6f  ",mymatrix.data[i]);
        if(mymatrix.data[i]!=0){
            printf("\n**test result of this function:FAIL**\n\n\n"); 
            break;
        }
        else{
            mymatrix_number++;
        }
        if(mymatrix_number==8){
             printf("\n**test result of this function:PASS**\n\n\n");
        }
    }

    printf("**testing is_null_matrix function**\n");
    printf("**param[in]- a pointer to a matrix**\n");
    printf("**param[out]-0(false) or 1(true)**\n");
    printf("**expected result:if it run successfully:is_null==0\n");
    printf("&&mymatrix_num==64&&mymatrix.data==0**\n");
    printf("**expected result:if it fails:is_null==1**\n");

    int is_null=is_null_matrix(&mymatrix);
    int mymatrix_num =get_matrix_size(&mymatrix);
    printf("**Factual result:is_null=%d\t",is_null);
    printf("mymatrix.data=%f\tmymatrix_num=%d**\n",*mymatrix.data,mymatrix_num);
    if(is_null==0){
         printf("**test result of this function:PASS**\n\n\n");
    }
    else{
         printf("**test result of this function:FAIL**\n\n\n");
    }
    
    d  = calculate_sd_matrix(sensor,n);

    double expect_d[8][8]={
    {1.000000,0.548810,0.606531,1.000000,0.670319,0.904839,0.030197,0.904835},
    {0.548810,1.000000,0.904835,0.548810,0.818730,0.496585,0.055023,0.606531},
    {0.606531,0.904835,1.000000,0.606531,0.904839,0.548812,0.049787,0.670322},
    {1.000000,0.548810,0.606531,1.000000,0.670319,0.904839,0.030197,0.904835},
    {0.670319,0.818730,0.904839,0.670319,1.000000,0.606531,0.045049,0.740819},
    {0.904839,0.496585,0.548812,0.904839,0.606531,1.000000,0.027324,0.818730},
    {0.030197,0.055023,0.049787,0.030197,0.045049,0.027324,1.000000,0.033373},
    {0.904835,0.606531,0.670322,0.904835,0.740819,0.818730,0.033373,1.000000}};
    int expect_value=0;
    printf("**testing calculate_sd_matrix function**\n");
    printf("**param[in]- an array of struct variable sensor and its numbers**\n");
    printf("**param[out]-the initial address of the struct variable sensor**\n");
    printf("**expected result:if it run successfully:expect_value==64 ");
    printf("and d[i][j]==expect_d[i][j]**\n");
    printf("**expected result:if it fails:have d[i][j]!=expect_d[i][j]**\n");
    printf("**expected value:");
    for(i=0;i<8;i++){
        printf("\n");
        for(j=0;j<8;j++){
            printf("%7.6f  ",d[i][j]);
            if(fabs(d[i][j]-expect_d[i][j])<0.0001){
                expect_value++;
            }
        }
    }
    printf("**");
    printf("\n**factual result:");
    for(i=0;i<8;i++){
        printf("\n");
        for(j=0;j<8;j++){
            printf("%7.6f  ",d[i][j]);
        }
    }
    printf("**\n");
    if(expect_value==64){
       printf("**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("**test result of this function:FAIL**\n\n\n");
    }
    for(i=0,k=0;i<n&&k<num;i++){
        for(j=0; j<n;j++){
            mymatrix.data[k] = d[i][j];
            k++;
        }
    }
    copy_matrix(&mymatrix, &temp);
    
    int test_copy=0;
    printf("**testing copy_matrix function**\n");
    printf("**param[in]- a pointer to a matrix being copied and a pointer"
        " to a matrix we get after copying**\n");
    printf("**param[out]-the initial address of the struct variable sensor**\n");
    printf("**expected result:if it run successfully :mymatrix.data[i] ==temp.data[i]"
        "&&mymatrix.column==temp.column&&mymatrix.row==temp.row\n");
    printf("**expected result:if it fails:have mymatrix.data[i]!=temp.data[i]"
           " or mymatrix.column!=temp.column or matrix.row!=temp.row**\n");
    printf("**Factual result:");
    printf("matrix being copied:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",temp.data[i]);
        if(mymatrix.data[i]==temp.data[i]&&mymatrix.column==
            temp.column&&mymatrix.row==temp.row){
            test_copy++;
        }
    }
    printf("\n**matrix we get after copying:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",mymatrix.data[i]);
    }
    if(test_copy==8){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }
    init_matrix(&temp_q, mymatrix.row, mymatrix.column);
    init_matrix(&temp_r, mymatrix.row, mymatrix.column);
    init_matrix(&eigen_value, mymatrix.row, 1);
    for (k=0;k<number_of_iterations;++k){
        qr_decomposition(&temp, &temp_q, &temp_r);
        matrix_mul_matrix(&temp, &temp_r, &temp_q);
    }
    for (k=0;k<temp.column;++k){
        eigen_value.data[k] = temp.data[k * temp.column + k];
    }
    
    sort_eigen_values(&eigen_value, 0);

    int eigen_values_compare=0;
    for(i=0;i<7;i++){
        if(eigen_value.data[i]>=eigen_value.data[i+1]){
            eigen_values_compare++;
        }
    }
    printf("**test sort_eigen_values function**\n");
    printf("**param[in]-a set of eigenvalues and a flag to decide whether"
        " in ascending or decending sequence**\n");
    printf("**param[out]-whether true(success) or false(fail)**\n");
    printf("**expected result:if it run successfully :eigen_values have been"
        " sorted in a descending sequence and the return value is 1**\n");
    printf("**expected result:if it fails: eigen_values have not been sorted"
        " in a descending sequence and the return value is 0**\n");
    printf("**Factual result:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",eigen_value.data[i]);
    }
    if(eigen_values_compare==7){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }
    calculate_eigenvector(&mymatrix, &temp_q, &eigen_value);
    
    eigenvectors_transpose_array = generate_eigenvector_transpose(&temp_q);
    init_matrix(&eigen_vectors, row, col);
    for(i=0,k=0;i<n&&k<num;i++){
        for(j=0; j<n;j++){
            eigen_vectors.data[k] = eigenvectors_transpose_array[i][j];
            k++;
        }
    }
    matrix_mul_matrix(&y, &eigen_vectors, &mymatrix);
    alpha_k = calculate_contribution_rate_of_kth_principal_component(&eigen_value,n);

    int check_alpha_value=0;
    double expect_alpha_value[8]=  {0.676500,0.138529,0.123636,0.023555,0.018072,0.011642,
        0.008066,0.000000};
    printf("**test sort_eigen_values function**\n");
    printf("**param[in]- the address of the first eigen_value and its number**\n");
    printf("**param[out]-the value of each principle component**\n");
    printf("**expected result:if it run successfully : the alpha_k is"
        " exactly what we expect**\n");
    printf("**expected result:if it fails: some of the principle component"
        " alpha_k is not what we expect **\n");
    printf("**Factual result:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",expect_alpha_value[i]);
    }
    printf("\n**alpha_k value:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",alpha_k[i]);
        if(fabs(expect_alpha_value[i]-alpha_k[i])<=1e-4){
            check_alpha_value++;
        }
    }
    if(check_alpha_value==8){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }

    for(i=0;i<n;i++){
        sum_varphi_m += alpha_k[i];
        m = i + 1;
        if(sum_varphi_m > p){
            break;
        }
    }

    z = compute_integrated_support_degree_score(&y,alpha_k,n,m);

    double expect_z_value[8]={1.513750,1.151085,1.242135,1.513750,1.312365,1.426167,
        0.160262,1.486718};
    int check_z_value=0;
    printf("**test comput_integrated_support_degree_score function**\n");
    printf("**param[in]- it is a pointer to a struct of matrix matrix\n"
           "*alpha_k - it is a pointer to an array of aplha_k\n"
           "n - it is a number of all sensors,\nm - it is a value that means we"
           " will select the first m principal components later**\n");
    printf("**param[out]-the support_degree_score z**\n");
    printf("**expected result:if it run successfully : the z is exactly what we expect\n");
    printf("**expected result:if it fails: some of the z value is not what we expect **\n");
    printf("**Factual result:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",z[i]);
    }
    printf("\n**z's actual value**:\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",expect_z_value[i]);
        if(fabs(expect_z_value[i]-z[i])<1e-4){
            check_z_value++;
        }
    }
    if(check_z_value==8){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }
    
    z_after = disregard_wrong_data(z,n,q);

    double expect_z_after[8]={1.513750,1.151085,1.242135,1.513750,1.312365,1.426167
        ,0.000000,1.486718
    };
    printf("**test comput_integrated_support_degree_score function **\n");
    printf("**param[in]- *z - it is a pointer to an array of all values of"
        " integrated support degree score\nn - it is a number of all sensors\nq"
        " - it is a changeable value that can limit the valid data**\n");
    printf("**param[out]-a pointer to an array of new integrated support degree"
        " score after disgarding wrong sensor data**\n");
    printf("**expected result:if it run successfully : the z_after is"
        " exactly what we expect\n");
    printf("**expected result:if it fails: some of the z_after value is"
        " not what we expect **\n");
    printf("**Factual result:**\n");
    int check_z_after=0;
    for(i=0;i<8;i++){
        printf("%7.6f   ",z_after[i]);
    }
    printf("\n**expect_z_after:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",expect_z_after[i]);
        if(fabs(expect_z_after[i]-z_after[i])<1e-4){
            check_z_after++;
        }
    }
    if(check_z_after==8){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }

    weight = compute_weight_coefficient(z_after,n);

    double expect_weight[8]={0.156931,0.119333,0.128772,0.156931,0.136053
        ,0.147851,0.000000,0.154128};
    int check_weight=0;
    printf("**test compute_weight_coefficient function**\n");
    printf("**param[in]- *z_after - it is a pointer to an array new integrated"
        " support degree score after disgarding wrong sensor data,n -"
        " it is a number of all sensors**\n");
    printf("**param[out]-a pointer to an array of weight coefficient of"
        " every valid sensor**\n");
    printf("**expected result:if it run successfully : the weight is exactly"
        " what we expected**\n");
    printf("**expected result:if it fails: some of the weight value is not"
        " what we expect **\n");
    printf("**Factual result:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",weight[i]);
    }
    printf("\n**expect_weight:**\n");
    for(i=0;i<8;i++){
        printf("%f   ",expect_weight[i]);
        if(fabs(weight[i]-expect_weight[i])<1e-4){
            check_weight++;
        }
    }
    if(check_weight==8){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }
    
    output = compute_fused_output(sensor,weight,n);

    double expected_output=53.008965;
    printf("**test compute_fused_output function**\n");
    printf("**param[in]- sensor - it is a struct of storing all sensor data,"
        "*weight - it is a pointer to an array of weight coefficient of every"
        " valid sensor,n - it is a number of all sensors**\n");
    printf("**param[out]-the value of fused output as the fianl result**\n");
    printf("**expected result:if it run successfully : the fused_output is"
        " exactly what we expected**\n");
    printf("**expected result:if it fails: the fused_output is"
        " not what we expected**\n");
    printf("**Factual result:%f**\n",expected_output);
    printf("**actual output:%f**\n",output);
    if(fabs(expected_output-output)<1e-4){
        printf("**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("**test result of this function:FAIL**\n\n\n");
    }
    
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

/**
 * \brief  to test the correctness of our whole algorithm and 
 *    each individual function component, and do some initialization and finishing touches
 * \details to manually compute each expected outputs of each function and compare it with 
 *    our actual outputs from the function and set some evaluation metrics to see 
 *    whether each of our function meets our standard
 * \param[in]  - the initial address of the struct variable sensor_list sensor_2
 * \return Success: true: each of function output meets our expected output and return a convincing result
 * \return Failure: some of the function output don't meet our expected output or 
 *    the process don't run normally ,or the expected final output is not what we get from the main program.
 */
double test_sensordataset2(struct sensor_list sensor[]){
    int i,j,k,m,n,row,col,num;
    n = 8;
    m = 0;
    row = n;
    col = n;
    num = row*col;
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
    matrix mymatrix,temp,temp_q,temp_r, eigen_value, eigen_vectors, y;
    init_matrix(&mymatrix, row, col);

    printf("**testing init_matrix function**\n");
    printf("**param[in]- a pointer to a matrix**\n");
    printf("**param[in]- the row and column of the matrix**\n");
    printf("**expected result:if it run successfully:row=8,column=8**\n");
    printf("**expected result:if it fails:row or column not equal to 8**\n");
    printf("**factual result:row=%d\tcolumn=%d**\n",mymatrix.row,mymatrix.column);
    if(mymatrix.column==8&&mymatrix.row==8){
        printf("**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("**test result of this function:FAIL**\n\n\n");
    }
    int size=get_matrix_size(&mymatrix);

    printf("**testing get_matrix_size function**\n");
    printf("**param[in]- a pointer to a matrix**\n");
    printf("**param[out]-size of the matrix**\n");
    printf("**expected result:if it run successfully:size=64**\n");
    printf("**expected result:if it fails:size not equal to 64**\n");
    printf("**Factual result: %d**\n",size);
    if(size==64){
        printf("**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("**test result of this function:FAIL**\n\n\n");
    }
    init_matrix(&temp, row, col);
    init_matrix(&y, row, col);
    set_matrix_zeros(&temp);
    set_matrix_zeros(&mymatrix);
    set_matrix_zeros(&y);

    int mymatrix_number=0;
    printf("**testing set_matrix zeros function**\n");
    printf("**param[in]- a pointer to a matrix**\n");
    printf("**param[out]-NULL**\n");
    printf("**expected result:if it run successfully:mymatrix.data[i]==0**\n");
    printf("**expected result:if it fails:have mymatrix.data[i]!=0**\n");
    printf("**Factual result:\n");
    for(i=0;i<8;i++){
        printf("%7.6f  ",mymatrix.data[i]);
        if(mymatrix.data[i]!=0) {
            printf("\n**test result of this function:FAIL**\n\n\n"); 
            break;
        }
        else{
            mymatrix_number++;
        }
        if(mymatrix_number==8){
             printf("\n**test result of this function:PASS**\n\n\n");
        }
    }

    printf("**testing is_null_matrix function**\n");
    printf("**param[in]- a pointer to a matrix**\n");
    printf("**param[out]-0(false) or 1(true)**\n");
    printf("**expected result:if it run successfully:is_null==0\n");
    printf("&&mymatrix_num==64&&mymatrix.data==0**\n");
    printf("**expected result:if it fails:is_null==1**\n");
    int is_null=is_null_matrix(&mymatrix);
    int mymatrix_num =get_matrix_size(&mymatrix);
    printf("**Factual result:is_null=%d\t",is_null);
    printf("mymatrix.data=%f\tmymatrix_num=%d**\n",*mymatrix.data,mymatrix_num);
    if(is_null==0){
         printf("**test result of this function:PASS**\n\n\n");
    }
    else{
         printf("**test result of this function:FAIL**\n\n\n");
    }
    
    d  = calculate_sd_matrix(sensor,n);

    double expect_d[8][8]={
    {1.000000,0.904837,0.670320,0.818731,0.818731,0.548812,0.067206,0.606531},
    {0.904837,1.000000,0.606531,0.904837,0.740818,0.606531,0.074274,0.670320},
    {0.670320,0.606531,1.000000,0.548812,0.818731,0.367879,0.045049,0.406570},
    {0.818731,0.904837,0.548812,1.000000,0.670320,0.670320,0.082085,0.740818},
    {0.818731,0.740818,0.818731,0.670320,1.000000,0.449329,0.055023,0.496585},
    {0.548812,0.606531,0.367879,0.670320,0.449329,1.000000,0.122456,0.904837},
    {0.067206,0.074274,0.045049,0.082085,0.055023,0.122456,1.000000,0.110803},
    {0.606531,0.670320,0.406570,0.740818,0.496585,0.904837,0.110803,1.000000}};
    int expect_value=0;
    printf("**testing calculate_sd_matrix function**\n");
    printf("**param[in]- an array of struct variable sensor and its numbers**\n");
    printf("**param[out]-the initial address of the struct variable sensor**\n");
    printf("**expected result:if it run successfully:expect_value==64 ");
    printf("and d[i][j]==expect_d[i][j]**\n");
    printf("**expected result:if it fails:have d[i][j]!=expect_d[i][j]**\n");
    printf("**expected value:");
    for(i=0;i<8;i++){
        printf("\n");
        for(j=0;j<8;j++){
            printf("%7.6f  ",d[i][j]);
            if(fabs(d[i][j]-expect_d[i][j])<0.0001){
                expect_value++;
            }
        }
    }
    printf("**");
    printf("\n**factual result:");
    for(i=0;i<8;i++){
        printf("\n");
        for(j=0;j<8;j++){
            printf("%7.6f  ",d[i][j]);
        }
    }
    printf("**\n");
    if(expect_value==64){
        printf("**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("**test result of this function:FAIL**\n\n\n");
    }
    for(i=0,k=0;i<n&&k<num;i++){
        for(j=0; j<n;j++){
            mymatrix.data[k] = d[i][j];
            k++;
        }
    }
    copy_matrix(&mymatrix, &temp);
    
    int test_copy=0;
    printf("**testing copy_matrix function**\n");
    printf("**param[in]- a pointer to a matrix being copied and a pointer"
        " to a matrix we get after copying**\n");
    printf("**param[out]-the initial address of the struct variable sensor**\n");
    printf("**expected result:if it run successfully :mymatrix.data[i] ==temp.data[i]"
        "&&mymatrix.column==temp.column&&mymatrix.row==temp.row\n");
    printf("**expected result:if it fails:have mymatrix.data[i]!=temp.data[i]"
           "or mymatrix.column!=temp.column or matrix.row!=temp.row**\n");
    printf("**Factual result:");
    printf("matrix being copied:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",temp.data[i]);
        if(mymatrix.data[i]==temp.data[i]&&mymatrix.column==temp.column&&mymatrix.row==temp.row){
            test_copy++;
        }
    }
    printf("\n**matrix we get after copying:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",mymatrix.data[i]);
    }
    if(test_copy==8){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }
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

    int eigen_values_compare=0;
    for(i=0;i<7;i++){
        if(eigen_value.data[i]>=eigen_value.data[i+1]){
            eigen_values_compare++;
        }
    }
    printf("**test sort_eigen_values function**\n");
    printf("**param[in]-a set of eigenvalues and a flag to decide whether"
        " in ascending or decending sequence**\n");
    printf("**param[out]-whether true(success) or false(fail)**\n");
    printf("**expected result:if it run successfully :eigen_values have been"
        " sorted in a descending sequence and the return value is 1**\n");
    printf("**expected result:if it fails: eigen_values have not been sorted"
        " in a descending sequence and the return value is 0**\n");
    printf("**Factual result:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",eigen_value.data[i]);
    }
    if(eigen_values_compare==7){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }
    
    
    calculate_eigenvector(&mymatrix, &temp_q, &eigen_value);

    eigenvectors_transpose_array = generate_eigenvector_transpose(&temp_q);
    init_matrix(&eigen_vectors, row, col);
    for(i=0,k=0;i<n&&k<num;i++){
        for(j=0; j<n;j++){
            eigen_vectors.data[k] = eigenvectors_transpose_array[i][j];
            k++;
        }
    }

    matrix_mul_matrix(&y, &eigen_vectors, &mymatrix);


    alpha_k = calculate_contribution_rate_of_kth_principal_component(&eigen_value,n);

    int check_alpha_value=0;
    double expect_alpha_value[8]=  {0.629353, 0.137459, 0.118052,
        0.059237, 0.022284, 0.014800, 0.011006, 0.007810};
    printf("**test sort_eigen_values function**\n");
    printf("**param[in]- the address of the first eigen_value and its number**\n");
    printf("**param[out]-the value of each principle component**\n");
    printf("**expected result:if it run successfully : the alpha_k is exactly"
        " what we expect**\n");
    printf("**expected result:if it fails: some of the principle component"
        " alpha_k is not what we expect **\n");
    printf("**Factual result:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",expect_alpha_value[i]);
    }
    printf("\n**alpha_k value:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",alpha_k[i]);
        if(fabs(expect_alpha_value[i]-alpha_k[i])<=1e-4){
            check_alpha_value++;
        }
    }
    if(check_alpha_value==8){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }
    for(i=0;i<n;i++){
        sum_varphi_m += alpha_k[i];
        m = i + 1;
        if(sum_varphi_m > p){
            break;
        }
    }
    
    z = compute_integrated_support_degree_score(&y,alpha_k,n,m);

    double expect_z_value[8]={1.256115,1.300495,0.952249,1.307302,1.114786,1.178317,
    0.152172,1.239587};
    int check_z_value=0;
    printf("**test comput_integrated_support_degree_score function**\n");
    printf("**param[in]- it is a pointer to a struct of matrix matrix\n"
           "*alpha_k - it is a pointer to an array of aplha_k\n"
           "n - it is a number of all sensors,\nm - it is a value that means"
           " we will select the first m principal components later**\n");
    printf("**param[out]-the support_degree_score z**\n");
    printf("**expected result:if it run successfully : the z is exactly"
        " what we expect\n ");
    printf("**expected result:if it fails: some of the z value is not what"
        " we expect **\n");
    printf("**Factual result:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",z[i]);
    }
    printf("\n**z's actual value**:\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",expect_z_value[i]);
        if(fabs(expect_z_value[i]-z[i])<1e-4){
            check_z_value++;
        }
    }
    if(check_z_value==8){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }
    
    z_after = disregard_wrong_data(z,n,q);

    double expect_z_after[8]={1.256115,1.300495,0.952249,1.307302,1.114786,1.178317
        ,0.000000,1.239587};
    printf("**test comput_integrated_support_degree_score function**\n");
    printf("**param[in]- *z - it is a pointer to an array of all values of "
        "integrated support degree score\nn - it is a number of all sensors\nq -"
        " it is a changeable value that can limit the valid data**\n");
    printf("**param[out]-a pointer to an array of new integrated support degree"
        " score after disgarding wrong sensor data**\n");
    printf("**expected result:if it run successfully : the z_after is"
        " exactly what we expect\n ");
    printf("**expected result:if it fails: some of the z_after value is"
        " not what we expect **\n");
    printf("**Factual result:**\n");
    int check_z_after=0;
    for(i=0;i<8;i++){
        printf("%7.6f   ",z_after[i]);
    }
    printf("\n**expect_z_after:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",expect_z_after[i]);
        if(fabs(expect_z_after[i]-z_after[i])<1e-4){
            check_z_after++;
        }
    }
    if(check_z_after==8){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }
    weight = compute_weight_coefficient(z_after,n);

    double expect_weight[8]={0.150454,0.155769,0.114057,0.156585,0.133526
    ,0.141135,0.000000,0.148474};
    int check_weight=0;
    printf("**test compute_weight_coefficient function**\n");
    printf("**param[in]- *z_after - it is a pointer to an array new integrated"
        " support degree score after disgarding wrong sensor data,n - "
        "it is a number of all sensors**\n");
    printf("**param[out]-a pointer to an array of weight coefficient of"
        " every valid sensor**\n");
    printf("**expected result:if it run successfully : the weight is exactly"
        " what we expected**\n ");
    printf("**expected result:if it fails: some of the weight value is not"
        " what we expect **\n");
    printf("**Factual result:**\n");
    for(i=0;i<8;i++){
        printf("%7.6f   ",weight[i]);
    }
    printf("\n**expect_weight:**\n");
    for(i=0;i<8;i++){
        printf("%f   ",expect_weight[i]);
        if(fabs(weight[i]-expect_weight[i])<1e-4){
            check_weight++;
        }
    }
    if(check_weight==8){
        printf("\n**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("\n**test result of this function:FAIL**\n\n\n");
    }
    
    output = compute_fused_output(sensor,weight,n);

    double expected_output=53.133484;
    printf("**test compute_fused_output function**\n");
    printf("**param[in]- sensor - it is a struct of storing all sensor data"
        ",*weight - it is a pointer to an array of weight coefficient of"
        " every valid sensor,n - it is a number of all sensors**\n");
    printf("**param[out]-the value of fused output as the fianl result**\n");
    printf("**expected result:if it run successfully : the fused_output is"
        " exactly what we expected**\n");
    printf("**expected result:if it fails: the fused_output is"
        " not what we expected**\n");
    printf("**Factual result:%f**\n",expected_output);
    printf("**actual output:%f**\n",output);
    if(fabs(expected_output-output)<1e-4){
        printf("**test result of this function:PASS**\n\n\n");
    }
    else{
        printf("**test result of this function:FAIL**\n\n\n");
    }
    
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

int main(){
    int choose = 0;
    printf("Please choose which set of data is going to be tested:"
            " (input 1 or 2, then press enter)\n");
    fflush(stdin);
    scanf("%d" ,&choose);

    struct sensor_list *sensor=input_data();
    struct sensor_list *sensor_1=&sensor[0];
    struct sensor_list *sensor_2=&sensor[8];

    switch(choose){
        case 1: test_sensordataset1(sensor_1); break;
        case 2: test_sensordataset2(sensor_2); break;
        default: printf("error! you need input right number!\n");
    }
    return 0;
}
