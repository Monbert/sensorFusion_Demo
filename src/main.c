/**
 * \brief Main file from which different functions located in other files are
 *    called to perform different tasks of the algorithm.
 * \author kaixia
 * \version 1.0
 * \date 2019-12-15
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "../include/algorithm.h"
#include "../include/data.h"


/** \brief The main function which distributes various tasks to other functions
 *  \return Success: 0
 *  \return Failure: NULL
 */
int main(){
    struct sensor_list *sensor=input_data();
    struct sensor_list *sensor_1=&sensor[0];
    struct sensor_list *sensor_2=&sensor[8];
    
    double output_1 = repeat(sensor_1);
    double output_2 = repeat(sensor_2);
    
    output_data(output_1,output_2);
    printf("Program running successfully. Results has been written in /data/output/outputdata.xls\n");

    return 0;
}
