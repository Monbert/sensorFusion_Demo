## sensorFusion_Demo<br>
Project of sensor fusion for course_sysc5709F<br>

Organization: Carleton University

Authors: <br>
Kai Xia   @Monbert<br>
Hui Tang  @htang085<br>

## Brief description of this software<br>
This software can calculate a single and weighted value from your input data which include many sensors values.

## File organization<br>
+ /bin  
 + executable
+ /build 
 + algorithm.o
 + data.o  
 + main.o
+ /data
 + /input/inputdata.csv
 + /output/outputdata.xls
+ /doc
 + software_design.pdf
 + software_instructions.pdf
+ /include
 + algorithm.h
 + data.h
+ /src
 + algorithm.c
 + data.c
 + main.c
+ /test
+ makefile
+ README

## Instructions<br>
Step 1 : Download the whole project <br>
Step 2 : cd to this project location in the terminal<br>
Step 3 : compile all files -- terminal command "make all"<br>
Step 4 : run the exectable file -- terminal command "make run"<br>
Step 5 : clean all files compiled -- terminal command "make clean"<br>
