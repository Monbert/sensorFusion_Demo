/**
 * \brief Function Prototype Library to Implement Data_input and Data_output
 * \details this file implements the input and output sensor data and 
 *    puts data into the sensorfuion algorithm
 * \author kaixia
 * \version 1.0
 * \date 2019-12-13
 */

#ifndef __DATA_H__
#define __DATA_H__


/**
 * \brief Function to output data into an .xls file
 * \details write the final outputs fused_value to the output file
 * \param[in] output_1 - the final fused_value
 * \param[in] output_2 - the final fused_value
 * \return Success: NULL
 * \return Failure: NULL
 */
void output_data(double output_1,double output_2);

/**
 * \brief Function to read the data from the input .csv file
 * \details read the initial data from sensors and store them into a strcut variable
 * \parameter[in] - NULL
 * \return Success: the initial address of the struct variable sensor
 * \return Failure: NULL
 */
struct sensor_list *input_data(void);

/**
 * \brief Function to set input sensor data into algorithm flow
 * \details repeatign to set input sensor data into algorithm flow and get results
 * \parameter[in] sensor - it is a struct of storing all sensor data
 * \return Success: the value of output sensor data
 * \return Failure: NULL
 */
double repeat(struct sensor_list sensor[]);


#endif //__DATA_H__
